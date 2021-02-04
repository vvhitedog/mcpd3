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
#include <memory>

#include "constraint.h"
#include <primaldual/mcpd3.h>
#include <graph/cycle.h>
#include <iostream>

namespace mcpd3 {

class DualDecomposition {
public:
  DualDecomposition(int npartition, int nnode, int narc,
                    std::vector<int> partitions, std::vector<int> arcs,
                    std::vector<int> arc_capacities,
                    std::vector<int> terminal_capacities)
      : npartition_(npartition), nnode_(nnode), narc_(narc),
        partitions_(std::move(partitions)), arcs_(std::move(arcs)), arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)), min_cut_sub_graphs_(npartition_), primal_solution_(nnode_) {
    initializeDecomposition();
  }

  DualDecomposition(int npartition, std::vector<int> partitions,
                    MinCutGraph min_cut_graph)
      : DualDecomposition(npartition, min_cut_graph.nnode, min_cut_graph.narc,
                          std::move(partitions), std::move(min_cut_graph.arcs),
                          std::move(min_cut_graph.arc_capacities),
                          std::move(min_cut_graph.terminal_capacities)) {}

  void runPrimalSolutionDecodingStep() {
    for ( int i = 0; i < npartition_; ++i ) {
      const auto &min_cut_sub_graph = min_cut_sub_graphs_[i];
      const auto &solver = solvers_[i];
      for ( const auto & [global_index,local_index] : min_cut_sub_graph.global_to_local_map ) {
        primal_solution_[global_index] = solver.getMinCutSolution(local_index);
      }
    }
    for ( auto &[global_index,constraints] : constraint_arc_map_ ) {
      int sum_x = 0;
      bool is_first_constraint = true;
      for ( auto &constraint : constraints ) {
        if ( is_first_constraint ) {
          sum_x += solvers_[constraint.partition_index_target].getMinCutSolution(constraint.local_index_target);
          is_first_constraint = false;
        }
        sum_x += solvers_[constraint.partition_index_source].getMinCutSolution(constraint.local_index_source);
      }
      primal_solution_[global_index]  = sum_x / (constraints.size() + 1); // use an averaging scheme to resolve what primal solution should be for constrained nodes
    }
  }


  long getPrimalMinCutValue() const {
      primal_solver_->setMinCutSolution(primal_solution_);
      return primal_solver_->getMinCutValue();
  }

  void runOptimizationStep(int nstep, int step_size, int max_cycle_count = 2) {
    CycleCountingList dual_cycle_list;
    long max_lower_bound = 0;
    long min_upper_bound = std::numeric_limits<long>::max();
    for (int i = 0; i < nstep; ++i) {

      long lower_bound = 0;
      long min_sub_lower_bound = std::numeric_limits<long>::max();
      long max_sub_lower_bound = std::numeric_limits<long>::min();
      for (auto &solver : solvers_) {
        solver.solve();
        auto adjusted_min_cut_value = solver.getMinCutValue();
        lower_bound += adjusted_min_cut_value;
        min_sub_lower_bound = std::min(min_sub_lower_bound,adjusted_min_cut_value);
        max_sub_lower_bound = std::max(max_sub_lower_bound,adjusted_min_cut_value);
        //printf(" -subproblem lower_bound- : %ld\n",adjusted_min_cut_value);
      }
      max_lower_bound = std::max(lower_bound,max_lower_bound);

      runPrimalSolutionDecodingStep();
      primal_solver_->setMinCutSolution(primal_solution_);
      auto upper_bound = primal_solver_->getMinCutValue();
      min_upper_bound = std::min(upper_bound,min_upper_bound);

      auto [num_disagreeing, disagreeing_global_indices] =
          runLagrangeMultipliersUpdateStep(step_size);
      printf("lower_bound : %ld num_disagreeing : %ld upper_bound : %ld\n",lower_bound,num_disagreeing,upper_bound);
      if (num_disagreeing == 0) { // optimality condition
        std::cout << "breaking because of no disagreement\n";
        break;
      }
      //std::stringstream stream_to_hash;
      //stream_to_hash << lower_bound << ";";
      //for (const auto &global_index : disagreeing_global_indices) {
      //  stream_to_hash << global_index << ";";
      //}
      //long hash_value = std::hash<std::string>{}(stream_to_hash.str());
      dual_cycle_list.addNode(getDualSolutionHash(std::move(disagreeing_global_indices),lower_bound));
      //dual_cycle_list.addNode(getDualSolutionHash(std::move(disagreeing_global_indices),0));
      if (dual_cycle_list.getMaxCycleCount() >
          max_cycle_count ) { // at least one set of configurations repeated twice
        std::cout << "breaking because cycle detected\n";
        break;
      }
  }
    printf(" === MAX === lower_bound : %ld\n",max_lower_bound);
    printf(" === MIN === upper_bound : %ld\n",min_upper_bound);
}

  template<int scale>
  void scaleProblem() {
      for (auto &solver : solvers_) {
        solver.scaleProblem<scale>();
      }
    for ( auto &[global_index,constraints] : constraint_arc_map_ ) {
      for ( auto &constraint : constraints ) {
        constraint.alpha *= scale;
      }
    }
  }


private:
  long getDualSolutionHash(std::list<int> disagreeing_global_indices, 
      long lower_bound ) const {
    std::hash<long> hasher{};
    long hash = hasher(lower_bound);
    for (const auto &global_index : disagreeing_global_indices ) {
      hash ^= hasher(global_index);
    }
    return hash;
  }

  std::pair<long,std::list<int>> runLagrangeMultipliersUpdateStep(int step_size) {
    std::list<int> disagreeing_global_indices;
    double num_disagreeing = 0;
    for ( auto &[global_index,constraints] : constraint_arc_map_ ) {
      for ( auto &constraint : constraints ) {
        int diff = solvers_[constraint.partition_index_target].getMinCutSolution(constraint.local_index_target)
          - solvers_[constraint.partition_index_source].getMinCutSolution(constraint.local_index_source);
        if ( diff != 0 ) {
        constraint.alpha += step_size*diff;
        disagreeing_global_indices.emplace_back(global_index);
        constraint.is_unsatisfied = 1;
        } else {
          constraint.is_unsatisfied = 0;
        }
        num_disagreeing += diff * diff;
      }
    }
    return {num_disagreeing,disagreeing_global_indices};
  }

  void initializeDecomposition() {
    std::unordered_map</*global_index=*/int,
                       /*exists_in_partitions=*/std::set<int>>
        constrained_nodes;
    /**
     * step 1: distribute all arcs into one and only one sub graph
     */
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      if ( s >  t ) {
        std::swap(s,t);
        std::swap(forward_capacity,backward_capacity);
      }
      // the sub graph each arc belongs to is defined to be the partition of the
      // source node
      int arc_partition = partitions_[s];
      auto &min_cut_sub_graph = min_cut_sub_graphs_[arc_partition];
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
      min_cut_sub_graphs_[partitions_[i]].insertTerminal(i, source_capacity,
                                                        sink_capacity);
    }
    /**
     * step 3: create solvers
     */
    for (auto &min_cut_sub_graph : min_cut_sub_graphs_) {
      solvers_.emplace_back(std::move(min_cut_sub_graph.graph));
    }
    for ( auto &solver : solvers_ )  {
      solver.setDecompositionSize(npartition_);
    }
    /**
     * step 4: create a DualDecompositionConstraintArc for each constraint
     * induced on each constrained node (which should be one less than the
     * number of partitions the node appears in)
     */
    std::set<int> constrained_nodes_partition_counts;
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
              min_cut_sub_graphs_[previous_partition].getNode(
                  global_index);
          int local_index_target =
              min_cut_sub_graphs_[partition].getNode(global_index);
          constraint_arcs.emplace_back(
              /*alpha=*/0,
              /*partition_index_source=*/previous_partition,
              /*partition_index_target=*/partition,
              /*local_index_source=*/local_index_source,
              /*local_index_target=*/local_index_target);
          auto arc_reference = --constraint_arcs.end();
          solvers_[previous_partition].addSourceDualDecompositionConstraint(
              arc_reference);
          solvers_[partition].addTargetDualDecompositionConstraint(arc_reference);
        }
        previous_partition = partition;
        not_first_partition = true;
      }
      assert(constraint_arcs.size() == partitions.size() - 1);
      constrained_nodes_partition_counts.insert(partitions.size());
    }
    printf("partition counts: ");
    for (const auto &count : constrained_nodes_partition_counts ) {
      printf("%d,",count);
    }
    printf("\n");
    /**
     * step 5: create a min cut problem from the original problem to evaluate the primal objective value
     */
    primal_solver_ = std::make_unique<PrimalDualMinCutSolver>(nnode_,narc_,
        std::move(arcs_),
        std::move(arc_capacities_),
        std::move(terminal_capacities_));
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

  template <typename T>
  /**
   * @brief a minimal class to make use of an expanding vector without the need
   * for copying its contents when expanding
   */
  class stable_vector {
  public:
    template <typename... Args> void emplace_back(Args &&... args) {
      list_.emplace_back(args...);
      vector_.emplace_back(--list_.end());
    }

    T &operator[](size_t index) { return *vector_[index]; }

    typename std::list<T>::iterator begin() {
      return list_.begin();
    }

    typename std::list<T>::iterator end() {
      return list_.end();
    }

  private:
    std::list<T> list_;
    std::vector<typename std::list<T>::iterator> vector_;
  };

  stable_vector<PrimalDualMinCutSolver> solvers_;
  std::unique_ptr<PrimalDualMinCutSolver> primal_solver_; // only used to evaluate primal value

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

    int getNode(int global_index) const {
      auto find_iter = global_to_local_map.find(global_index);
      if (find_iter == global_to_local_map.end()) {
        throw std::runtime_error("Node not found");
      }
      return find_iter->second;
    }

    void insertArc(int global_source_index, int global_target_index,
                   int forward_capacity, int backward_capacity) {
      int s = getOrInsertNode(global_source_index);
      int t = getOrInsertNode(global_target_index);
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

  std::vector<MinCutSubGraph> min_cut_sub_graphs_;
  std::vector<int> primal_solution_;
};

} // namespace mcpd3
