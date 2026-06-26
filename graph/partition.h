// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <atomic>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <limits>
#include <graph/mcgraph.h>
#include <string>
#include <thread>
#include <vector>
#ifdef HAVE_METIS
#include <metis.h>
#endif

namespace mcpd3 {

inline bool partition_progress_enabled() {
  const char *value = std::getenv("MCPD3_PROGRESS");
  return value != nullptr && value[0] != '\0' && std::string(value) != "0";
}

inline void partition_progress_report(
    const char *stage, long done, long total,
    std::chrono::steady_clock::time_point start) {
  if (!partition_progress_enabled()) {
    return;
  }
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
  const double pct = total > 0 ? 100.0 * static_cast<double>(done) / total : 0;
  const double rate = elapsed > 0 ? static_cast<double>(done) / elapsed : 0;
  const double eta = rate > 0 && total > done
                         ? static_cast<double>(total - done) / rate
                         : 0;
  std::fprintf(stderr,
               "mcpd3_progress stage=%s done=%ld total=%ld pct=%.2f "
               "elapsed_sec=%.1f eta_sec=%.1f rate=%.0f_per_sec\n",
               stage, done, total, pct, elapsed, eta, rate);
  std::fflush(stderr);
}

std::vector<int> basic_graph_partition(int npartition, int narc, int nnode,
                                       const std::vector<int> &arc) {
  std::vector<int> partitions(nnode);
  int nnode_in_each_partition = (nnode + npartition - 1) / npartition;
  for (int i = 0; i < nnode; ++i) {
    partitions[i] = i / nnode_in_each_partition;
  }
  return std::move(partitions);
}

std::vector<int> basic_graph_partition(int npartition,
                                       const mcpd3::MinCutGraph &graph) {
  return basic_graph_partition(npartition, graph.narc, graph.nnode, graph.arcs);
}

inline int env_int_or_default(const char *name, int default_value) {
  const char *value = std::getenv(name);
  return value != nullptr && value[0] != '\0' ? std::atoi(value)
                                              : default_value;
}

inline double env_double_or_default(const char *name, double default_value) {
  const char *value = std::getenv(name);
  return value != nullptr && value[0] != '\0' ? std::atof(value)
                                              : default_value;
}

inline double partition_edge_weight(const std::vector<int> *arc_capacities,
                                    int edge_index, double lambda) {
  if (arc_capacities == nullptr || arc_capacities->empty() || lambda == 0.0) {
    return 1.0;
  }
  const size_t base = static_cast<size_t>(2) * edge_index;
  const double cap0 = base < arc_capacities->size() ? (*arc_capacities)[base] : 0;
  const double cap1 =
      base + 1 < arc_capacities->size() ? (*arc_capacities)[base + 1] : 0;
  const double cap = std::max(cap0, cap1);
  return 1.0 + lambda * std::log1p(std::max(0.0, cap));
}

inline bool partition_is_boundary_node(int node, int part,
                                       const std::vector<int> &partitions,
                                       const std::vector<int> &xadj,
                                       const std::vector<int> &adjv) {
  for (int idx = xadj[node]; idx < xadj[node + 1]; ++idx) {
    if (partitions[adjv[idx]] != part) {
      return true;
    }
  }
  return false;
}

inline std::vector<int> region_grow_initial_partition(
    int npartition, int nnode, const std::vector<int> &xadj,
    const std::vector<int> &adjv, int max_part_size) {
  const bool report_progress = partition_progress_enabled();
  const long progress_interval = 10000000;
  const auto start = std::chrono::steady_clock::now();
  std::vector<int> partitions(nnode, -1);
  std::vector<int> part_sizes(npartition, 0);
  std::vector<std::deque<int>> queues(npartition);

  for (int p = 0; p < npartition; ++p) {
    const int seed =
        npartition == 1
            ? 0
            : static_cast<int>((static_cast<long long>(p) * (nnode - 1)) /
                               (npartition - 1));
    if (partitions[seed] == -1) {
      partitions[seed] = p;
      ++part_sizes[p];
      for (int idx = xadj[seed]; idx < xadj[seed + 1]; ++idx) {
        queues[p].push_back(adjv[idx]);
      }
    }
  }

  int next_unassigned = 0;
  long assigned = 0;
  for (int size : part_sizes) {
    assigned += size;
  }
  bool changed = true;
  while (assigned < nnode && changed) {
    changed = false;
    for (int p = 0; p < npartition && assigned < nnode; ++p) {
      if (part_sizes[p] >= max_part_size) {
        continue;
      }
      while (!queues[p].empty() && partitions[queues[p].front()] != -1) {
        queues[p].pop_front();
      }
      int node = -1;
      if (!queues[p].empty()) {
        node = queues[p].front();
        queues[p].pop_front();
      } else {
        while (next_unassigned < nnode && partitions[next_unassigned] != -1) {
          ++next_unassigned;
        }
        if (next_unassigned < nnode) {
          node = next_unassigned;
        }
      }
      if (node == -1 || partitions[node] != -1) {
        continue;
      }
      partitions[node] = p;
      ++part_sizes[p];
      ++assigned;
      changed = true;
      for (int idx = xadj[node]; idx < xadj[node + 1]; ++idx) {
        if (partitions[adjv[idx]] == -1) {
          queues[p].push_back(adjv[idx]);
        }
      }
      if (report_progress && assigned % progress_interval == 0) {
        partition_progress_report("local_partition_region_grow", assigned,
                                  nnode, start);
      }
    }
  }

  for (int node = 0; node < nnode; ++node) {
    if (partitions[node] != -1) {
      continue;
    }
    int best_part = 0;
    for (int p = 1; p < npartition; ++p) {
      if (part_sizes[p] < part_sizes[best_part]) {
        best_part = p;
      }
    }
    partitions[node] = best_part;
    ++part_sizes[best_part];
  }
  partition_progress_report("local_partition_region_grow", nnode, nnode,
                            start);
  return partitions;
}

inline std::vector<int> local_search_graph_partition(
    int npartition, int narc, int nnode, const std::vector<int> &arc,
    const std::vector<int> *arc_capacities = nullptr) {
  const bool report_progress = partition_progress_enabled();
  const long progress_interval = 10000000;
  const int passes = std::max(0, env_int_or_default("MCPD3_LOCAL_PARTITION_PASSES", 3));
  const double lambda = env_double_or_default("MCPD3_LOCAL_PARTITION_LAMBDA", 0.005);
  const double balance_slack =
      std::max(0.0, env_double_or_default("MCPD3_LOCAL_PARTITION_BALANCE_SLACK", 0.05));
  const auto partition_start = std::chrono::steady_clock::now();

  std::vector<int> xadj(nnode + 1);
  auto degree_start = std::chrono::steady_clock::now();
  for (int aid = 0; aid < narc; ++aid) {
    const int s = arc[2 * aid + 0];
    const int t = arc[2 * aid + 1];
    ++xadj[s + 1];
    ++xadj[t + 1];
    if (report_progress && (aid + 1) % progress_interval == 0) {
      partition_progress_report("local_partition_degree_count", aid + 1, narc,
                                degree_start);
    }
  }
  partition_progress_report("local_partition_degree_count", narc, narc,
                            degree_start);

  auto prefix_start = std::chrono::steady_clock::now();
  for (int i = 0; i < nnode; ++i) {
    xadj[i + 1] += xadj[i];
    if (report_progress && (i + 1) % progress_interval == 0) {
      partition_progress_report("local_partition_prefix_sum", i + 1, nnode,
                                prefix_start);
    }
  }
  partition_progress_report("local_partition_prefix_sum", nnode, nnode,
                            prefix_start);

  std::vector<int> adjv(xadj[nnode]);
  std::vector<float> adjw(xadj[nnode]);
  std::vector<int> cursor = xadj;
  auto fill_start = std::chrono::steady_clock::now();
  for (int aid = 0; aid < narc; ++aid) {
    const int s = arc[2 * aid + 0];
    const int t = arc[2 * aid + 1];
    const double weight = partition_edge_weight(arc_capacities, aid, lambda);
    int pos = cursor[s]++;
    adjv[pos] = t;
    adjw[pos] = weight;
    pos = cursor[t]++;
    adjv[pos] = s;
    adjw[pos] = weight;
    if (report_progress && (aid + 1) % progress_interval == 0) {
      partition_progress_report("local_partition_adjacency_fill", aid + 1,
                                narc, fill_start);
    }
  }
  partition_progress_report("local_partition_adjacency_fill", narc, narc,
                            fill_start);
  cursor.clear();
  cursor.shrink_to_fit();

  const int target_size = (nnode + npartition - 1) / npartition;
  const int max_part_size =
      std::max(1, static_cast<int>(std::ceil(target_size * (1.0 + balance_slack))));
  const char *init_env = std::getenv("MCPD3_LOCAL_PARTITION_INIT");
  const std::string init_mode = init_env == nullptr ? "region" : init_env;
  std::vector<int> partitions =
      init_mode == "contiguous" || init_mode == "basic"
          ? basic_graph_partition(npartition, narc, nnode, arc)
          : region_grow_initial_partition(npartition, nnode, xadj, adjv,
                                          max_part_size);
  std::vector<int> part_sizes(npartition, 0);
  for (int p : partitions) {
    ++part_sizes[p];
  }

  if (report_progress) {
    std::fprintf(stderr,
                 "mcpd3_progress stage=local_partition_start nnode=%d "
                 "narc=%d partitions=%d passes=%d init=%s lambda=%.6g "
                 "balance_slack=%.6g max_part_size=%d\n",
                 nnode, narc, npartition, passes, init_mode.c_str(), lambda,
                 balance_slack, max_part_size);
    std::fflush(stderr);
  }

  std::vector<int> neighbor_parts;
  std::vector<double> neighbor_weights;
  for (int pass = 0; pass < passes; ++pass) {
    auto pass_start = std::chrono::steady_clock::now();
    long moves = 0;
    double improvement = 0.0;
    for (int node = 0; node < nnode; ++node) {
      const int old_part = partitions[node];
      neighbor_parts.clear();
      neighbor_weights.clear();

      double old_internal_weight = 0.0;
      for (int idx = xadj[node]; idx < xadj[node + 1]; ++idx) {
        const int nbr_part = partitions[adjv[idx]];
        const double weight = adjw[idx];
        if (nbr_part == old_part) {
          old_internal_weight += weight;
          continue;
        }
        auto it = std::find(neighbor_parts.begin(), neighbor_parts.end(),
                            nbr_part);
        if (it == neighbor_parts.end()) {
          neighbor_parts.push_back(nbr_part);
          neighbor_weights.push_back(weight);
        } else {
          neighbor_weights[static_cast<size_t>(it - neighbor_parts.begin())] +=
              weight;
        }
      }

      int best_part = old_part;
      double best_delta = 0.0;
      const bool old_node_boundary = partition_is_boundary_node(
          node, old_part, partitions, xadj, adjv);
      for (size_t cand = 0; cand < neighbor_parts.size(); ++cand) {
        const int new_part = neighbor_parts[cand];
        if (new_part == old_part || part_sizes[new_part] >= max_part_size) {
          continue;
        }

        double delta = old_internal_weight - neighbor_weights[cand];

        const bool new_node_boundary = partition_is_boundary_node(
            node, new_part, partitions, xadj, adjv);
        if (old_node_boundary && !new_node_boundary) {
          delta -= 1.0;
        } else if (!old_node_boundary && new_node_boundary) {
          delta += 1.0;
        }

        for (int idx = xadj[node]; idx < xadj[node + 1]; ++idx) {
          const int nbr = adjv[idx];
          const int nbr_part = partitions[nbr];
          const bool before_boundary = partition_is_boundary_node(
              nbr, nbr_part, partitions, xadj, adjv);
          partitions[node] = new_part;
          const bool after_boundary = partition_is_boundary_node(
              nbr, nbr_part, partitions, xadj, adjv);
          partitions[node] = old_part;
          if (before_boundary && !after_boundary) {
            delta -= 1.0;
          } else if (!before_boundary && after_boundary) {
            delta += 1.0;
          }
        }

        if (delta < best_delta) {
          best_delta = delta;
          best_part = new_part;
        }
      }

      if (best_part != old_part) {
        partitions[node] = best_part;
        --part_sizes[old_part];
        ++part_sizes[best_part];
        ++moves;
        improvement -= best_delta;
      }

      if (report_progress && (node + 1) % progress_interval == 0) {
        partition_progress_report("local_partition_pass", node + 1, nnode,
                                  pass_start);
      }
    }
    partition_progress_report("local_partition_pass", nnode, nnode,
                              pass_start);
    if (report_progress) {
      std::fprintf(stderr,
                   "mcpd3_progress stage=local_partition_pass_summary "
                   "pass=%d moves=%ld improvement=%.6f\n",
                   pass + 1, moves, improvement);
      std::fflush(stderr);
    }
    if (moves == 0) {
      break;
    }
  }

  if (report_progress) {
    long boundary_nodes = 0;
    double crossing_weight = 0.0;
    for (int node = 0; node < nnode; ++node) {
      if (partition_is_boundary_node(node, partitions[node], partitions, xadj,
                                     adjv)) {
        ++boundary_nodes;
      }
    }
    for (int aid = 0; aid < narc; ++aid) {
      const int s = arc[2 * aid + 0];
      const int t = arc[2 * aid + 1];
      if (partitions[s] != partitions[t]) {
        crossing_weight += partition_edge_weight(arc_capacities, aid, lambda);
      }
    }
    std::fprintf(stderr,
                 "mcpd3_progress stage=local_partition_summary "
                 "boundary_nodes=%ld crossing_weight=%.6f\n",
                 boundary_nodes, crossing_weight);
    std::fflush(stderr);
  }
  partition_progress_report("local_partition_total", 1, 1, partition_start);

  return partitions;
}

#ifdef HAVE_METIS

std::vector<int> metis_partition(int npartition, int narc, int nnode,
                                 const std::vector<int> &arc) {

  const bool report_progress = partition_progress_enabled();
  const long progress_interval = 10000000;
  const auto partition_start = std::chrono::steady_clock::now();

  std::vector<idx_t> xadj(nnode + 1);
  auto degree_start = std::chrono::steady_clock::now();
  for (size_t aid = 0; aid < static_cast<size_t>(narc); ++aid) {
    const int32_t s = arc[2 * aid + 0];
    const int32_t t = arc[2 * aid + 1];
    ++xadj[s + 1];
    ++xadj[t + 1];
    if (report_progress && (aid + 1) % progress_interval == 0) {
      partition_progress_report("metis_degree_count", aid + 1, narc,
                                degree_start);
    }
  }
  partition_progress_report("metis_degree_count", narc, narc, degree_start);
  auto prefix_start = std::chrono::steady_clock::now();
  for (int inode = 0; inode < nnode; ++inode) {
    xadj[inode + 1] += xadj[inode];
    if (report_progress && (inode + 1) % progress_interval == 0) {
      partition_progress_report("metis_prefix_sum", inode + 1, nnode,
                                prefix_start);
    }
  }
  partition_progress_report("metis_prefix_sum", nnode, nnode, prefix_start);

  std::vector<idx_t> adjv(static_cast<size_t>(xadj[nnode]));
  std::vector<idx_t> cursor = xadj;
  auto fill_start = std::chrono::steady_clock::now();
  for (size_t aid = 0; aid < static_cast<size_t>(narc); ++aid) {
    const int32_t s = arc[2 * aid + 0];
    const int32_t t = arc[2 * aid + 1];
    adjv[cursor[s]++] = t;
    adjv[cursor[t]++] = s;
    if (report_progress && (aid + 1) % progress_interval == 0) {
      partition_progress_report("metis_adjacency_fill", aid + 1, narc,
                                fill_start);
    }
  }
  partition_progress_report("metis_adjacency_fill", narc, narc, fill_start);
  cursor.clear();
  cursor.shrink_to_fit();

  std::vector<idx_t> part(nnode);
  idx_t ncon = 1;
  idx_t _nnode = nnode;
  idx_t _npart = npartition;
  idx_t objval;
  int ret = METIS_ERROR;
  auto metis_start = std::chrono::steady_clock::now();
  if (report_progress) {
    std::fprintf(stderr,
                 "mcpd3_progress stage=metis_call_start nnode=%d narc=%d "
                 "adjacency_entries=%ld partitions=%d\n",
                 nnode, narc, static_cast<long>(xadj[nnode]), npartition);
    std::fflush(stderr);
  }
  if (report_progress) {
    std::atomic<bool> metis_done(false);
    std::thread metis_thread([&] {
      ret = METIS_PartGraphKway(&_nnode, &ncon, &xadj[0], &adjv[0], nullptr,
                                nullptr, nullptr, &_npart, nullptr, nullptr,
                                nullptr, &objval, &part[0]);
      metis_done.store(true, std::memory_order_release);
    });
    while (!metis_done.load(std::memory_order_acquire)) {
      std::this_thread::sleep_for(std::chrono::seconds(5));
      const double elapsed =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        metis_start)
              .count();
      if (!metis_done.load(std::memory_order_acquire)) {
        std::fprintf(stderr,
                     "mcpd3_progress stage=metis_call elapsed_sec=%.1f "
                     "eta_sec=unknown reason=opaque_metis_call\n",
                     elapsed);
        std::fflush(stderr);
      }
    }
    metis_thread.join();
  } else {
    ret = METIS_PartGraphKway(&_nnode, &ncon, &xadj[0], &adjv[0], nullptr,
                              nullptr, nullptr, &_npart, nullptr, nullptr,
                              nullptr, &objval, &part[0]);
  }
  partition_progress_report("metis_call_done", 1, 1, metis_start);
  partition_progress_report("metis_partition_total", 1, 1, partition_start);

  if (ret != METIS_OK) {
    throw std::runtime_error("Something failed while partitioning.\n");
  }

  std::vector<int> label(nnode, -1);
  std::copy(part.begin(), part.end(), label.begin());
  return std::move(label);
}

std::vector<int> metis_partition(int npartition,
                                 const mcpd3::MinCutGraph &graph) {
  const auto &narc = graph.narc;
  const auto &arc = graph.arcs;
  const auto &nnode = graph.nnode;
  return metis_partition(npartition, narc, nnode, arc);
}
#endif

inline std::vector<int> configured_graph_partition(
    int npartition, int narc, int nnode, const std::vector<int> &arc,
    const std::vector<int> *arc_capacities = nullptr) {
  const char *mode_env = std::getenv("MCPD3_PARTITIONER");
  const std::string mode = mode_env == nullptr ? "metis" : mode_env;
  if (mode == "basic" || mode == "contiguous") {
    return basic_graph_partition(npartition, narc, nnode, arc);
  }
  if (mode == "local" || mode == "local_search") {
    return local_search_graph_partition(npartition, narc, nnode, arc,
                                        arc_capacities);
  }
#ifdef HAVE_METIS
  if (mode != "metis") {
    std::fprintf(stderr,
                 "mcpd3_progress stage=partitioner_warning unknown=%s "
                 "fallback=metis\n",
                 mode.c_str());
    std::fflush(stderr);
  }
  return metis_partition(npartition, narc, nnode, arc);
#else
  if (mode != "basic" && mode != "contiguous") {
    std::fprintf(stderr,
                 "mcpd3_progress stage=partitioner_warning unknown=%s "
                 "fallback=basic\n",
                 mode.c_str());
    std::fflush(stderr);
  }
  return basic_graph_partition(npartition, narc, nnode, arc);
#endif
}
} // namespace mcpd3
