/* graph.cpp */

#include "graph.h"
#include <errno.h>
#include <fcntl.h>
#include <new>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <type_traits>
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

mcpd3::SolverStorageOptions get_bk_storage_options() {
  mcpd3::SolverStorageOptions options;
  const char *mode = getenv("MCPD3_BK_STORAGE");
  if (mode && mode[0] != '\0') {
    if (strcmp(mode, "malloc") == 0) {
      options.mode = mcpd3::SolverStorageMode::RESIDENT;
      return options;
    }
    if (strcmp(mode, "file_mmap") == 0) {
      options.mode = mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
      const char *directory = getenv("MCPD3_BK_MMAP_DIR");
      options.directory = directory ? directory : "";
      const char *advice = getenv("MCPD3_BK_MMAP_ADVISE");
      options.mmap_advice = advice ? advice : "";
      return options;
    }
    if (strcmp(mode, "anon_mmap") == 0 || strcmp(mode, "anonymous_mmap") == 0) {
      options.mode = mcpd3::SolverStorageMode::ANONYMOUS_MMAP;
      const char *advice = getenv("MCPD3_BK_MMAP_ADVISE");
      options.mmap_advice = advice ? advice : "";
      return options;
    }
    fprintf(stderr, "unknown MCPD3_BK_STORAGE=%s; using malloc\n", mode);
    return options;
  }

  const char *mmap_dir = getenv("MCPD3_BK_MMAP_DIR");
  if (mmap_dir && mmap_dir[0] != '\0') {
    options.mode = mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
    options.directory = mmap_dir;
  }
  const char *advice = getenv("MCPD3_BK_MMAP_ADVISE");
  options.mmap_advice = advice ? advice : "";
  return options;
}

void *allocate_bk_array(size_t bytes, const char *kind, int &fd,
                        bool &is_mmap_backed, bool &is_file_backed,
                        const mcpd3::SolverStorageOptions &storage_options) {
  fd = -1;
  is_mmap_backed = false;
  is_file_backed = false;
  if (storage_options.mode == mcpd3::SolverStorageMode::RESIDENT) {
    return malloc(bytes);
  }

  int mmap_flags = MAP_SHARED;
  if (storage_options.mode == mcpd3::SolverStorageMode::ANONYMOUS_MMAP) {
    mmap_flags = MAP_PRIVATE | MAP_ANONYMOUS;
#ifdef MAP_POPULATE
    if (storage_options.mmap_advice == "populate") {
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
    try {
      mcpd3::applySolverMmapAdvice(ptr, bytes, storage_options, kind);
    } catch (...) {
      munmap(ptr, bytes);
      throw;
    }
    return ptr;
  }

  if (storage_options.directory.empty()) {
    fprintf(stderr,
            "file-backed BK storage requires an explicit directory\n");
    return nullptr;
  }

  std::string pattern =
      storage_options.directory + "/mcpd3_bk_" + kind + "_XXXXXX";
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
  is_file_backed = true;
  try {
    mcpd3::applySolverMmapAdvice(ptr, bytes, storage_options, kind);
  } catch (...) {
    munmap(ptr, bytes);
    close(fd);
    fd = -1;
    throw;
  }
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
    : Graph(node_num_max, edge_num_max, get_bk_storage_options(), err_function) {}

template <typename captype, typename tcaptype, typename flowtype>
Graph<captype, tcaptype, flowtype>::Graph(
    int node_num_max, int edge_num_max,
    const mcpd3::SolverStorageOptions &storage_options,
    void (*err_function)(const char *))
    : nodes_mmap_backed(false), arcs_mmap_backed(false),
      nodes_file_backed(false), arcs_file_backed(false), nodes_mmap_fd(-1),
      arcs_mmap_fd(-1), nodes_mmap_bytes(0), arcs_mmap_bytes(0), node_num(0),
      nodeptr_block(NULL), error_function(err_function) {
  if (node_num_max < 16)
    node_num_max = 16;
  if (edge_num_max < 16)
    edge_num_max = 16;

  if (storage_options.mode != mcpd3::SolverStorageMode::RESIDENT) {
    changed_arc_marks = mcpd3::SolverArray<unsigned char>(
        static_cast<std::size_t>(edge_num_max), static_cast<unsigned char>(0),
        storage_options, "bk_changed_arc_marks");
  }

  nodes_mmap_bytes = node_num_max * sizeof(node);
  arcs_mmap_bytes = 2 * edge_num_max * sizeof(arc);
  if constexpr (std::is_trivially_copyable_v<node>) {
    nodes = (node *)allocate_bk_array(nodes_mmap_bytes, "nodes", nodes_mmap_fd,
                                      nodes_mmap_backed, nodes_file_backed,
                                      storage_options);
  } else {
    nodes = new (std::nothrow) node[static_cast<size_t>(node_num_max)];
  }
  if constexpr (std::is_trivially_copyable_v<arc>) {
    arcs = (arc *)allocate_bk_array(arcs_mmap_bytes, "arcs", arcs_mmap_fd,
                                    arcs_mmap_backed, arcs_file_backed,
                                    storage_options);
  } else {
    arcs = new (std::nothrow) arc[static_cast<size_t>(2 * edge_num_max)];
  }
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
  if constexpr (std::is_trivially_copyable_v<node>) {
    free_bk_array(nodes, nodes_mmap_bytes, nodes_mmap_fd, nodes_mmap_backed);
  } else {
    delete[] nodes;
  }
  if constexpr (std::is_trivially_copyable_v<arc>) {
    free_bk_array(arcs, arcs_mmap_bytes, arcs_mmap_fd, arcs_mmap_backed);
  } else {
    delete[] arcs;
  }
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
  if constexpr (std::is_trivially_copyable_v<node>) {
    nodes = (node *)realloc(nodes_old, node_num_max * sizeof(node));
  } else {
    nodes = new (std::nothrow) node[static_cast<size_t>(node_num_max)];
    if (nodes) {
      for (int index = 0; index < node_num; ++index) {
        nodes[index] = nodes_old[index];
      }
    }
  }
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
  if constexpr (!std::is_trivially_copyable_v<node>) {
    delete[] nodes_old;
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
  if constexpr (std::is_trivially_copyable_v<arc>) {
    arcs = (arc *)realloc(arcs_old, arc_num_max * sizeof(arc));
  } else {
    arcs = new (std::nothrow) arc[static_cast<size_t>(arc_num_max)];
    if (arcs) {
      for (int index = 0; index < arc_num; ++index) {
        arcs[index] = arcs_old[index];
      }
    }
  }
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
  if constexpr (!std::is_trivially_copyable_v<arc>) {
    delete[] arcs_old;
  }
}

#include "instances.inc"
