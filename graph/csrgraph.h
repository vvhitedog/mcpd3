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

#include <cassert>
#include <string>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include <graph/dimacs.h>
#include <graph/partition.h>
#include <io/mmaparray.h>
#include <io/workdir.h>
#include <maxflow/graph.h>

namespace mcpd3 {

template <typename node_index_type, typename arc_index_type, typename cap_type>
struct Edge {
  node_index_type source, target;
  cap_type capacity;
  bool operator<(const Edge &rhs) const {
    return source < rhs.source || (source == rhs.source && target < rhs.target);
  }

  bool operator==(const Edge &rhs) const {
    return source == rhs.source && target == rhs.target;
  }
};

template <typename node_index_type = int, typename arc_index_type = long,
          typename cap_type = int, typename flow_type = long>
class CsrGraph {
public:
  CsrGraph(const std::string file_dir_prefix)
      : file_dir_prefix_(file_dir_prefix) {
    create_work_directory(file_dir_prefix_);
  }

  void setCut( const std::vector<bool> &cut ) {
    cut_ = cut;
  }

  flow_type getCurrentCutValue() const {
    flow_type cut_value = 0;
    for ( node_index_type inode = 0; inode < nnode_; ++inode ) {
      // 1. calculate arc contribution
      arc_index_type arc_start_offset = (*adjacency_offsets_)[inode];
      arc_index_type arc_end_offset = (*adjacency_offsets_)[inode+1];
      for ( arc_index_type iarc = arc_start_offset; iarc < arc_end_offset; ++iarc ) {
        auto end_node = (*adjacency_nodes_)[iarc];
        auto arc_cap = (*arc_capacities_)[iarc];
        if ( !cut_[inode] && cut_[end_node] ) {
          cut_value += arc_cap;
        }
      }
      
      // 2. calculate terminal contribution
      auto terminal_cap = (*terminal_capacities_)[inode];
      if ( terminal_cap > 0 && cut_[inode] ) {
        cut_value += terminal_cap;
      }
      if ( terminal_cap < 0 && !cut_[inode] ) {
        cut_value += -terminal_cap;
      }
    }
    return cut_value;
  }

  void narrowBandDecode(const std::list<node_index_type> &seeds,
      node_index_type max_distance = 14 ){

    // 1. run bfs to get nodes within proximity of seeds
    std::unordered_map<node_index_type,node_index_type> distance_map;
    std::unordered_set<node_index_type> visited;
    std::queue<node_index_type> to_visit;
    arc_index_type arcs_visited = 0;

    for ( const auto &seed : seeds ) {
      distance_map[seed] = 0;
      visited.insert(seed);
      to_visit.push(seed);
    }

    while ( !to_visit.empty() ){
      auto current_node = to_visit.front();
      to_visit.pop();
      auto current_dist = distance_map[current_node];
      if ( current_dist >= max_distance ) {
        continue;
      }
      arc_index_type arc_start_offset = (*adjacency_offsets_)[current_node];
      arc_index_type arc_end_offset = (*adjacency_offsets_)[current_node+1];
      for ( arc_index_type iarc = arc_start_offset; iarc < arc_end_offset; ++iarc, ++arcs_visited ) {
        auto end_node = (*adjacency_nodes_)[iarc];
        auto [_,success] = visited.insert(end_node);
        if ( success ) { // already visited
          to_visit.push(end_node);
          distance_map[end_node] = current_dist + 1;
        }
      }
    }

    // 2. remap indices
    std::unordered_map<node_index_type,node_index_type> new_index;
    node_index_type index = 0;
    for ( const auto &node : visited ) {
      new_index[node] = index++;
    }
    node_index_type node_count = visited.size();

    // 3. create a subset maxflow graph

    // 3 a. gather arcs
    std::unordered_map<std::pair<node_index_type,node_index_type>,
      std::pair<cap_type,cap_type>,boost::hash<std::pair<node_index_type,node_index_type>> >
      arc_to_capacity_map;
    for ( const auto &start_node: visited ) {
      arc_index_type arc_start_offset = (*adjacency_offsets_)[start_node];
      arc_index_type arc_end_offset = (*adjacency_offsets_)[start_node+1];
      for ( arc_index_type iarc = arc_start_offset; iarc < arc_end_offset; ++iarc ) {
        auto end_node = (*adjacency_nodes_)[iarc];
        if ( visited.find(end_node) == visited.end() ) { // skip this arc
          continue;
        }
        auto arc_cap = (*arc_capacities_)[iarc];
        if ( start_node > end_node ) {
          auto key = std::make_pair(end_node,start_node);
          arc_to_capacity_map[key].second = arc_cap;
        } else {
          auto key = std::make_pair(start_node,end_node);
          arc_to_capacity_map[key].first = arc_cap;
        }
      }
    }

    // 3 b. create maxflow graph and add arcs
    MaxflowGraph maxflow_graph(node_count,arc_to_capacity_map.size());
    maxflow_graph.add_node(node_count);
    for ( const auto &[st,caps] : arc_to_capacity_map ) {
      const auto &[source,target] = st;
      const auto &[forward_cap,backward_cap] = caps;
      maxflow_graph.add_edge(new_index[source],new_index[target],forward_cap,backward_cap);
    }
    arc_to_capacity_map.clear();

    // 3 c. create terminals nodes
    for ( const auto &[node,distance] : distance_map ) {
      if ( distance < max_distance ) {
        auto terminal_cap = (*terminal_capacities_)[node];
        maxflow_graph.set_trcap(new_index[node],terminal_cap);
      } else {
        if ( !cut_[node] ) {
          maxflow_graph.set_trcap(new_index[node],std::numeric_limits<cap_type>::max()/2);
        }
        if ( cut_[node] ) {
          maxflow_graph.set_trcap(new_index[node],-std::numeric_limits<cap_type>::max()/2);
        }
      }
    }

    // 4. solve maxflow
    maxflow_graph.maxflow();

    // 5. decode solution for visited nodes
    for ( const auto &node : visited ) {
      cut_[node] = maxflow_graph.what_segment(new_index[node]) == MaxflowGraph::SINK;
    }

  }

  template <typename SortedEdgesArray>
  void initializeArcsFromSortedEdges(node_index_type nnode, arc_index_type narc,
                                     const SortedEdgesArray &sorted_edges) {
    nnode_ = nnode;
    narc_ = narc;
    adjacency_offsets_ = std::make_unique<MmapArray<arc_index_type>>(
        nnode_ + 1, file_dir_prefix_ + "/adjacency_offsets.mmap");
    adjacency_nodes_ = std::make_unique<MmapArray<node_index_type>>(
        narc_, file_dir_prefix_ + "/adjacency_nodes.mmap");
    arc_capacities_ = std::make_unique<MmapArray<cap_type>>(
        narc_, file_dir_prefix_ + "/arc_capacities.mmap");

    node_index_type current_node = -1;
    arc_index_type arc_processed_count = 0;
    for (auto edge_it = sorted_edges.begin(); arc_processed_count < narc_;
         ++edge_it) {
      if (edge_it->source != current_node) {
        for ( auto i = current_node; i < edge_it->source; ++i ) {
          (*adjacency_offsets_)[i + 1] = arc_processed_count;
        }
        current_node = edge_it->source;
        if (current_node >= nnode_) {
          throw std::runtime_error("current node is larger than expected");
        }
      }
      (*arc_capacities_)[arc_processed_count] = edge_it->capacity;
      (*adjacency_nodes_)[arc_processed_count++] = edge_it->target;
    }
    (*adjacency_offsets_)[nnode_] = arc_processed_count;
  }

  void initializeTerminalCapacities(
      std::unique_ptr<MmapArray<cap_type>> terminal_capacities) {
    terminal_capacities_ = std::move(terminal_capacities);
  }

#ifdef HAVE_METIS
  void partitionWithMetis(node_index_type npartition) {
    {
    MmapArray<idx_t> adjacency_offsets(nnode_+1,file_dir_prefix_ + "/adjacency_offsets_metis.mmap");
    MmapArray<idx_t> adjacency_nodes(narc_,file_dir_prefix_ + "/adjacency_nodes_metis.mmap");
    std::copy(adjacency_offsets_->begin(),adjacency_offsets_->end(),adjacency_offsets.begin());
    std::copy(adjacency_nodes_->begin(),adjacency_nodes_->end(),adjacency_nodes.begin());
    MmapArray<idx_t> partitions(nnode_,file_dir_prefix_ + "/partitions_metis.mmap");

    idx_t ncon = 1;
    idx_t _nnode = nnode_;
    idx_t _npart = npartition;
    idx_t objval;
    int ret = METIS_PartGraphKway(&_nnode, &ncon, adjacency_offsets.data(), adjacency_nodes.data(), nullptr,
        nullptr, nullptr, &_npart, nullptr, nullptr,
        nullptr, &objval, partitions.data());

    parition_labels_ = std::make_unique<MmapArray<cap_type>>(
        nnode_, file_dir_prefix_ + "/partition_labels.mmap");
    std::copy(partitions.begin(),partitions.end(),parition_labels_->begin());
    }

    // clean up unused memmaps
    std::filesystem::remove(file_dir_prefix_ + "/adjacency_offsets_metis.mmap");
    std::filesystem::remove(file_dir_prefix_ + "/adjacency_nodes_metis.mmap");
    std::filesystem::remove(file_dir_prefix_ + "/partitions_metis.mmap");
  }
#endif

private:
  node_index_type nnode_;
  arc_index_type narc_;
  std::string file_dir_prefix_;

  std::unique_ptr<MmapArray<arc_index_type>> adjacency_offsets_;
  std::unique_ptr<MmapArray<node_index_type>> adjacency_nodes_;
  std::unique_ptr<MmapArray<cap_type>> arc_capacities_;
  std::unique_ptr<MmapArray<cap_type>> terminal_capacities_;
  std::unique_ptr<MmapArray<node_index_type>> parition_labels_;

  std::vector<bool> cut_; // stored in memory due to small memory footprint


  using MaxflowGraph =
    Graph</*captype=*/cap_type, /*tcaptype=*/cap_type, /*flowtype=*/flow_type>;
};

template <typename node_index_type = int, typename arc_index_type = long,
          typename cap_type = int>
CsrGraph<node_index_type, arc_index_type, cap_type>
read_dimacs_to_csr(const std::string &filename,
                   const std::string &work_dir = "./csr_graph_dimacs") {
  // 0. create work directory
  create_work_directory(work_dir);

  // 1. read once through to get edge (upper-bound) and node count
  arc_index_type edge_count = 0;
  node_index_type num_nodes = 0;
  {
    auto arc_op = [&](node_index_type s, node_index_type t, cap_type cap) {
      num_nodes = std::max<int>(num_nodes, s + 1);
      num_nodes = std::max<int>(num_nodes, t + 1);
      edge_count += 2;
    };
    auto term_op = [&](bool is_source, node_index_type n, cap_type cap) {
      num_nodes = std::max<int>(num_nodes, n + 1);
    };
    _dimacs_implementation::read_dimacs_general(filename, arc_op, term_op);
  }

  // 2. read in the edges & terminal capacities (to be used later)
  using Edge = Edge<node_index_type, arc_index_type, cap_type>;
  MmapArray<Edge> edges(edge_count, work_dir + "/edges.mmap");
  std::unique_ptr<MmapArray<cap_type>> terminal_capacities =
      std::make_unique<MmapArray<cap_type>>(
          num_nodes, work_dir + "/terminal_capacities.mmap");
  std::fill(terminal_capacities->begin(), terminal_capacities->end(), 0);
  edge_count = 0;
  {
    auto arc_op = [&](node_index_type s, node_index_type t, cap_type cap) {
      edges[edge_count++] = Edge{s, t, cap};
      edges[edge_count++] = Edge{t, s, 0};
    };
    auto term_op = [&](bool is_source, node_index_type n, cap_type cap) {
      (*terminal_capacities)[n] += is_source ? cap : -cap;
    };
    _dimacs_implementation::read_dimacs_general(filename, arc_op, term_op);
  }

  // 3. sort the edges
  std::sort(edges.begin(), edges.end());

  // 4. accumulate capacities of duplicate edges and remove duplicates
  Edge *current_edge = nullptr;
  for (auto &edge : edges) {
    if (current_edge == nullptr || !(edge == *current_edge)) {
      current_edge = &edge;
    } else {
      current_edge->capacity += edge.capacity;
      edge.source = std::numeric_limits<node_index_type>::max();
      edge.target = 0;
    }
  }
  std::sort(edges.begin(), edges.end());
  auto end_arcs =
      std::lower_bound(edges.begin(), edges.end(),
                       Edge{std::numeric_limits<node_index_type>::max(), 0, 0});
  auto num_edges = std::distance(edges.begin(), end_arcs);

  // 5. create the CSR graph from the sorted edges and terminal capacities
  CsrGraph<node_index_type, arc_index_type, cap_type> g(work_dir);
  g.initializeArcsFromSortedEdges(num_nodes, num_edges, edges);
  g.initializeTerminalCapacities(std::move(terminal_capacities));

  return std::move(g);
}

} // namespace mcpd3
