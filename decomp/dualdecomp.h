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

#include <list>
#include <set>
#include <unordered_map>

#include "constraint.h"
#include <primaldual/mcpd3.h>

namespace mcpd3 {

class DualDecomposition {
public:
  DualDecomposition(int npartition, int nnode, int narc,
                    std::vector<int> &&partitions, std::vector<int> &&arcs,
                    std::vector<int> &&arc_capacities,
                    std::vector<int> &&terminal_capacities)
      : npartition_(npartition), nnode_(nnode), narc_(narc),
        partitions_(partitions), arcs_(arcs), arc_capacities_(arc_capacities),
        terminal_capacities_(terminal_capacities) {
          initializeDecomposition();
        }

private:
  void initializeDecomposition() {
    std::unordered_map</*global_index=*/int,
                       /*exists_in_partitions=*/std::set<int>>
        constrained_nodes;
    std::vector<MinCutSubGraph> min_cut_sub_graphs(npartition_);
    /**
     * step 1: distribute all arcs into one and only one sub graph
     */
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      // the sub graph each arc belongs to is defined to be the partition of the
      // source node
      int arc_partition = partitions_[s];
      auto &min_cut_sub_graph = min_cut_sub_graphs[arc_partition];
      min_cut_sub_graph.insertArc(s, t, forward_capacity, backward_capacity);
      if (partitions_[t] != arc_partition) { // t is an auxillary node that
                                             // needs to be constrained
        constrained_nodes[t].insert(arc_partition);
      }
    }
    /**
     * step 2: add source and sink capacities of nodes
     */
    for (int i = 0; i < nnode_; ++i) {
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      min_cut_sub_graphs[partitions_[i]].insertTerminal(i, source_capacity,
                                                        sink_capacity);
    }
    /**
     * step 3: create solvers
     */
    for (auto &min_cut_sub_graph : min_cut_sub_graphs) {
      solvers_.emplace_back(min_cut_sub_graph.graph);
    }
    /**
     * step 4: create DualDecompositionConstraintArc for each constraint
     * induced on each constrained node (which should be one less than the
     * number of partitions the node appears in)
     */
    for (auto &[global_index, partitions] : constrained_nodes) {
      partitions.insert(
          partitions_[global_index]); // list each constrained node in its
                                      // original partition
      auto &constraint_arcs = constraint_arc_map_[global_index];
      assert(partitions.size() >
             1); // requirement to be a proper constrained node
      bool not_first_partition = false;
      int previous_partition = -1;
      for (const auto &partition :
           partitions) { // iterates in sorted order due to std::set
        if (not_first_partition) {
          int local_index_source =
              min_cut_sub_graphs[previous_partition].getOrInsertNode(
                  global_index);
          int local_index_target =
              min_cut_sub_graphs[partition].getOrInsertNode(global_index);
          constraint_arcs.emplace_back(
              /*alpha=*/0,
              /*partition_index_source=*/previous_partition,
              /*partition_index_target=*/partition,
              /*local_index_source=*/local_index_source,
              /*local_index_target=*/local_index_target);
          auto arc_reference = --constraint_arcs.end();
          solvers_[previous_partition].addDualDecompositionConstraint(
              arc_reference, true);
          solvers_[partition].addDualDecompositionConstraint(arc_reference,
                                                             false);
        }
        previous_partition = partition;
        not_first_partition = true;
      }
      assert(constraint_arcs.size() == partitions.size() - 1);
    }
  }

  /**
   * data passed into decomposition
   */
  int npartition_;
  int nnode_;
  int narc_;
  std::vector<int> partitions_;
  std::vector<int> arcs_;
  std::vector<int> arc_capacities_;
  std::vector<int> terminal_capacities_;

  /**
   * data structures needed for solving dual decomposition
   */
  std::unordered_map</*global_index=*/int,
                     std::list<DualDecompositionConstraintArc>>
      constraint_arc_map_;
  std::vector<PrimalDualMinCutSolver> solvers_;

  struct MinCutSubGraph {
    MinCutGraph graph;
    std::unordered_map</*global_index=*/int, /*local_index=*/int>
        global_to_local_map;

    MinCutSubGraph() {
      graph.nnode = 0;
      graph.narc = 0;
    }

    int getOrInsertNode(int global_index) {
      auto find_iter = global_to_local_map.find(global_index);
      if (find_iter == global_to_local_map.end()) {
        auto [insert_iter, exists] =
            global_to_local_map.insert({global_index, graph.nnode++});
        return insert_iter->second;
      }
      return find_iter->second;
    }

    void insertArc(int global_source_index, int global_target_index,
                   int forward_capacity, int backward_capacity) {
      int s = getOrInsertNode(global_source_index);
      int t = getOrInsertNode(global_target_index);
      assert(s < t);
      graph.arc_capacities.push_back(forward_capacity);
      graph.arc_capacities.push_back(backward_capacity);
      graph.arcs.push_back(s);
      graph.arcs.push_back(t);
      graph.narc++;
    }

    void insertTerminal(int global_index, int source_capacity,
                        int sink_capacity) {
      int s = getOrInsertNode(global_index);
      graph.terminal_capacities.resize(2 * graph.nnode, 0);
      graph.terminal_capacities[2 * s + 0] = source_capacity;
      graph.terminal_capacities[2 * s + 1] = sink_capacity;
    }
  };
};

} // namespace mcpd3
