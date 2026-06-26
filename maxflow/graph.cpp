/* graph.cpp */

#include "graph.h"
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/*
        special constants for node->parent. Duplicated in maxflow.cpp, both
   should match!
*/
#define TERMINAL ((arc *)1) /* to terminal */
#define ORPHAN ((arc *)2)   /* orphan */

namespace {

enum class BKStorageMode {
  Malloc,
  FileMmap,
  AnonymousMmap,
};

BKStorageMode get_bk_storage_mode() {
  const char *mode = getenv("MCPD3_BK_STORAGE");
  if (mode && mode[0] != '\0') {
    if (strcmp(mode, "malloc") == 0) {
      return BKStorageMode::Malloc;
    }
    if (strcmp(mode, "file_mmap") == 0) {
      return BKStorageMode::FileMmap;
    }
    if (strcmp(mode, "anon_mmap") == 0 || strcmp(mode, "anonymous_mmap") == 0) {
      return BKStorageMode::AnonymousMmap;
    }
    fprintf(stderr, "unknown MCPD3_BK_STORAGE=%s; using malloc\n", mode);
    return BKStorageMode::Malloc;
  }

  const char *mmap_dir = getenv("MCPD3_BK_MMAP_DIR");
  if (mmap_dir && mmap_dir[0] != '\0') {
    return BKStorageMode::FileMmap;
  }
  return BKStorageMode::Malloc;
}

void advise_bk_mapping(void *ptr, size_t bytes, const char *kind) {
  const char *advise = getenv("MCPD3_BK_MMAP_ADVISE");
  if (!advise || advise[0] == '\0' || strcmp(advise, "none") == 0) {
    return;
  }
  int rc = 0;
  if (strcmp(advise, "willneed") == 0) {
    rc = madvise(ptr, bytes, MADV_WILLNEED);
#ifdef MAP_POPULATE
  } else if (strcmp(advise, "populate") == 0) {
    // MAP_POPULATE is applied at mmap time. Keep this spelling accepted so
    // users can combine one environment interface with both mmap sites.
    rc = madvise(ptr, bytes, MADV_WILLNEED);
#endif
#ifdef MLOCK_ONFAULT
  } else if (strcmp(advise, "lock_onfault") == 0) {
    rc = mlock2(ptr, bytes, MLOCK_ONFAULT);
#endif
  } else if (strcmp(advise, "lock") == 0) {
    rc = mlock(ptr, bytes);
  } else if (strcmp(advise, "dontdump") == 0) {
#ifdef MADV_DONTDUMP
    rc = madvise(ptr, bytes, MADV_DONTDUMP);
#endif
  } else {
    fprintf(stderr, "unknown MCPD3_BK_MMAP_ADVISE=%s; ignoring\n", advise);
    return;
  }
  if (rc != 0) {
    fprintf(stderr, "BK mmap advise %s failed for %s array of %zu bytes: %s\n",
            advise, kind, bytes, strerror(errno));
  }
}

void *allocate_bk_array(size_t bytes, const char *kind, int &fd,
                        bool &is_mmap_backed) {
  fd = -1;
  is_mmap_backed = false;
  const BKStorageMode mode = get_bk_storage_mode();
  if (mode == BKStorageMode::Malloc) {
    return malloc(bytes);
  }

  int mmap_flags = MAP_SHARED;
  if (mode == BKStorageMode::AnonymousMmap) {
    mmap_flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef MAP_POPULATE
    const char *advise = getenv("MCPD3_BK_MMAP_ADVISE");
    if (advise && strcmp(advise, "populate") == 0) {
      mmap_flags |= MAP_POPULATE;
    }
#endif
    void *ptr = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, mmap_flags, -1, 0);
    if (ptr == MAP_FAILED) {
      fprintf(stderr, "failed to anonymous mmap BK %s array of %zu bytes: %s\n",
              kind, bytes, strerror(errno));
      return nullptr;
    }
    is_mmap_backed = true;
    advise_bk_mapping(ptr, bytes, kind);
    return ptr;
  }

  const char *mmap_dir = getenv("MCPD3_BK_MMAP_DIR");
  if (!mmap_dir || mmap_dir[0] == '\0') {
    fprintf(stderr,
            "MCPD3_BK_STORAGE=file_mmap requires MCPD3_BK_MMAP_DIR\n");
    return nullptr;
  }

  std::string pattern = std::string(mmap_dir) + "/mcpd3_bk_" + kind + "_XXXXXX";
  fd = mkstemp(pattern.data());
  if (fd == -1) {
    fprintf(stderr, "failed to create BK mmap file %s: %s\n", pattern.c_str(),
            strerror(errno));
    return nullptr;
  }
  unlink(pattern.c_str());
  if (ftruncate(fd, bytes) != 0) {
    fprintf(stderr, "failed to size BK mmap file %s to %zu bytes: %s\n",
            pattern.c_str(), bytes, strerror(errno));
    close(fd);
    fd = -1;
    return nullptr;
  }
  void *ptr = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (ptr == MAP_FAILED) {
    fprintf(stderr, "failed to mmap BK %s array of %zu bytes: %s\n", kind,
            bytes, strerror(errno));
    close(fd);
    fd = -1;
    return nullptr;
  }
  is_mmap_backed = true;
  advise_bk_mapping(ptr, bytes, kind);
  return ptr;
}

void free_bk_array(void *ptr, size_t bytes, int fd, bool is_mmap_backed) {
  if (is_mmap_backed) {
    if (ptr) {
      munmap(ptr, bytes);
    }
    if (fd != -1) {
      close(fd);
    }
  } else {
    free(ptr);
  }
}

} // namespace

template <typename captype, typename tcaptype, typename flowtype>
Graph<captype, tcaptype, flowtype>::Graph(int node_num_max, int edge_num_max,
                                          void (*err_function)(const char *))
    : nodes_mmap_backed(false), arcs_mmap_backed(false), nodes_mmap_fd(-1),
      arcs_mmap_fd(-1), nodes_mmap_bytes(0), arcs_mmap_bytes(0), node_num(0),
      nodeptr_block(NULL), error_function(err_function) {
  if (node_num_max < 16)
    node_num_max = 16;
  if (edge_num_max < 16)
    edge_num_max = 16;

  nodes_mmap_bytes = node_num_max * sizeof(node);
  arcs_mmap_bytes = 2 * edge_num_max * sizeof(arc);
  nodes = (node *)allocate_bk_array(nodes_mmap_bytes, "nodes", nodes_mmap_fd,
                                    nodes_mmap_backed);
  arcs = (arc *)allocate_bk_array(arcs_mmap_bytes, "arcs", arcs_mmap_fd,
                                  arcs_mmap_backed);
  if (!nodes || !arcs) {
    if (error_function)
      (*error_function)("Not enough memory!");
    exit(1);
  }

  node_last = nodes;
  node_max = nodes + node_num_max;
  arc_last = arcs;
  arc_max = arcs + 2 * edge_num_max;

  maxflow_iteration = 0;
  flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype>
Graph<captype, tcaptype, flowtype>::~Graph() {
  if (nodeptr_block) {
    delete nodeptr_block;
    nodeptr_block = NULL;
  }
  free_bk_array(nodes, nodes_mmap_bytes, nodes_mmap_fd, nodes_mmap_backed);
  free_bk_array(arcs, arcs_mmap_bytes, arcs_mmap_fd, arcs_mmap_backed);
}

template <typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reset() {
  node_last = nodes;
  arc_last = arcs;
  node_num = 0;

  if (nodeptr_block) {
    delete nodeptr_block;
    nodeptr_block = NULL;
  }

  maxflow_iteration = 0;
  flow = 0;
}

template <typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reallocate_nodes(int num) {
  if (nodes_mmap_backed) {
    if (error_function)
      (*error_function)("BK mmap mode does not support node reallocation");
    exit(1);
  }
  int node_num_max = (int)(node_max - nodes);
  node *nodes_old = nodes;

  node_num_max += node_num_max / 2;
  if (node_num_max < node_num + num)
    node_num_max = node_num + num;
  nodes = (node *)realloc(nodes_old, node_num_max * sizeof(node));
  if (!nodes) {
    if (error_function)
      (*error_function)("Not enough memory!");
    exit(1);
  }

  node_last = nodes + node_num;
  node_max = nodes + node_num_max;

  if (nodes != nodes_old) {
    node *i;
    arc *a;
    for (i = nodes; i < node_last; i++) {
      if (i->next)
        i->next =
            (node *)((char *)i->next + (((char *)nodes) - ((char *)nodes_old)));
    }
    for (a = arcs; a < arc_last; a++) {
      a->head =
          (node *)((char *)a->head + (((char *)nodes) - ((char *)nodes_old)));
    }
  }
}

template <typename captype, typename tcaptype, typename flowtype>
void Graph<captype, tcaptype, flowtype>::reallocate_arcs() {
  if (arcs_mmap_backed) {
    if (error_function)
      (*error_function)("BK mmap mode does not support arc reallocation");
    exit(1);
  }
  int arc_num_max = (int)(arc_max - arcs);
  int arc_num = (int)(arc_last - arcs);
  arc *arcs_old = arcs;

  arc_num_max += arc_num_max / 2;
  if (arc_num_max & 1)
    arc_num_max++;
  arcs = (arc *)realloc(arcs_old, arc_num_max * sizeof(arc));
  if (!arcs) {
    if (error_function)
      (*error_function)("Not enough memory!");
    exit(1);
  }

  arc_last = arcs + arc_num;
  arc_max = arcs + arc_num_max;

  if (arcs != arcs_old) {
    node *i;
    arc *a;
    for (i = nodes; i < node_last; i++) {
      if (i->first)
        i->first =
            (arc *)((char *)i->first + (((char *)arcs) - ((char *)arcs_old)));
      if (i->parent && i->parent != ORPHAN && i->parent != TERMINAL)
        i->parent =
            (arc *)((char *)i->parent + (((char *)arcs) - ((char *)arcs_old)));
    }
    for (a = arcs; a < arc_last; a++) {
      if (a->next)
        a->next =
            (arc *)((char *)a->next + (((char *)arcs) - ((char *)arcs_old)));
      a->sister =
          (arc *)((char *)a->sister + (((char *)arcs) - ((char *)arcs_old)));
    }
  }
}

#include "instances.inc"
