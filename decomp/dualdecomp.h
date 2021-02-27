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

#include <cmath>
#include <list>
#include <memory>
#include <queue>
#include <set>
#include <unordered_map>
#include <iostream>

#include <decomp/constraint.h>
#include <multithread/threadpool.h>
#include <graph/cycle.h>
#include <graph/partition.h>
#include <primaldual/mcpd3.h>

namespace mcpd3 {

/**
 * @brief Time a lambda function.
 *
 * @param lambda - the function to execute and time
 *
 * @return the number of microseconds elapsed while executing lambda
 */
template <typename Lambda>
std::chrono::microseconds time_lambda(Lambda lambda) {
  auto start_time = std::chrono::high_resolution_clock::now();
  lambda();
  auto end_time = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                               start_time);
}

class DualDecomposition {
public:
  DualDecomposition(int npartition, int nnode, int narc,
                    std::vector<int> arcs,
                    std::vector<int> arc_capacities,
                    std::vector<int> terminal_capacities)
      : npartition_(npartition), nnode_(nnode), narc_(narc),
        arcs_(std::move(arcs)),
        arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)),
        min_cut_sub_graphs_(npartition_), primal_solution_(nnode_), scale_(1),
        thread_pool_(std::min<size_t>(npartition_,std::thread::hardware_concurrency()) ), solve_loop_time_(0),
        max_lower_bound_(std::numeric_limits<double>::lowest()) {
    initializeDecomposition();
  }

  DualDecomposition(int npartition,
                    MinCutGraph min_cut_graph)
      : DualDecomposition(npartition, min_cut_graph.nnode, min_cut_graph.narc,
                          std::move(min_cut_graph.arcs),
                          std::move(min_cut_graph.arc_capacities),
                          std::move(min_cut_graph.terminal_capacities)) {}

  long getTotalSolveLoopTime() const { return solve_loop_time_; }

  void runPrimalSolutionDecodingStep(bool do_narrow_band_decode = false) {
    for (int i = 0; i < npartition_; ++i) {
      const auto &min_cut_sub_graph = min_cut_sub_graphs_[i];
      const auto &solver = solvers_[i];
      for (const auto &[global_index, local_index] :
           min_cut_sub_graph.global_to_local_map) {
        primal_solution_[global_index] = solver->getMinCutSolution(local_index);
      }
    }
    int total_disagree_count = 0;
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      double sum_x = 0;
      bool disagreement = false;
      for (auto &constraint : constraints) {
        auto u = solvers_[constraint.partition_index_target]->getMinCutSolution(
            constraint.local_index_target);
        auto v = solvers_[constraint.partition_index_source]->getMinCutSolution(
            constraint.local_index_source);
        if (u != v) {
          disagreement = true;
        }
        sum_x += u;
        sum_x += v;
      }
      primal_solution_[global_index] = std::round(
          sum_x /
          (2 * constraints
                   .size())); // use an averaging scheme to resolve what primal
                              // solution should be for constrained nodes
      if (disagreement) {
        total_disagree_count++;
      }
    }

    if (do_narrow_band_decode && total_disagree_count > 0) {
      std::list<int> disagree_nodes;
      for (auto &[global_index, constraints] : constraint_arc_map_) {
        bool disagreement = false;
        for (auto &constraint : constraints) {
          auto u =
              solvers_[constraint.partition_index_target]->getMinCutSolution(
                  constraint.local_index_target);
          auto v =
              solvers_[constraint.partition_index_source]->getMinCutSolution(
                  constraint.local_index_source);
          if (u != v) {
            disagreement = true;
            break;
          }
        }
        if (disagreement) { // search for better configuration
          disagree_nodes.emplace_back(global_index);
        }
      }
      primal_solver_->setMinCutSolution(primal_solution_);
      primal_solver_->decodeNarrowBand(disagree_nodes, 14);
      printf("recalc primal: %ld \n", primal_solver_->getMinCutValue());
    }
  }

  enum OptimizationStatus {
    OPTIMAL,
    NO_FURTHER_PROGRESS,
    ITERATION_COUNT_EXCEEDED
  };

  template<bool attempt_decoding, typename Decoder>
  void solve(Decoder decoder) {
    const int scaling_factor = 10;
    const int num_optimization_scales = 5;
    const int max_cycle_count = 2;
    for (int iscale = 0; iscale < num_optimization_scales; ++iscale) {
      OptimizationStatus status;
      auto ms = time_lambda([&] {
      status = runOptimizationScale(10000, 1, max_cycle_count, true);
      });
      std::cout << " opt scale runtime ms: " << ms.count() << "\n";
      if (status == mcpd3::DualDecomposition::OPTIMAL) {
        break;
      }
      if ( attempt_decoding ) { // decoding requested
        if ( decoder(primal_solution_,max_lower_bound_,disagreeing_global_indices_) ) {
          break;
        }
      }
      auto ms2 = time_lambda([&] {
      scaleProblem<scaling_factor>();
      });
      std::cout << " scale prob runtime ms: " << ms2.count() << "\n";
    }
  }

  void solve() {
    auto null_decoder = [=](const std::vector<bool> &cut,
        double max_lower_bound,
        const std::list<int> &disagreeing_global_indices ) -> bool{
        return false;
        };
    solve<false>(null_decoder);
  }

  OptimizationStatus runOptimizationScale(int nstep, int step_size,
                                          int max_cycle_count = 2,
                                          int break_on_small_change = false) {
    OptimizationStatus opt_status = ITERATION_COUNT_EXCEEDED;
    const int num_stats_in_group = 10;
    TwoGroupScalarStatisticsTracker<long> lower_bound_group_stats(
        num_stats_in_group);
    CycleCountingList dual_cycle_list;
    long max_lower_bound = 0;
    for (int i = 0; i < nstep; ++i) {

      std::atomic<long> lower_bound = 0;
      auto solve_loop_time = time_lambda([&] {
        for (auto &solver : solvers_) {
          thread_pool_.push([&] {
            solver->solve();
            lower_bound += solver->getMinCutValue();
          });
        }
        thread_pool_.wait();
      });
      solve_loop_time_ += solve_loop_time.count();
      max_lower_bound =
          std::max(static_cast<long>(lower_bound), max_lower_bound);

      disagreeing_global_indices_ =
          runLagrangeMultipliersUpdateStep(step_size, break_on_small_change);
      printf("lower_bound : %lf num_disagreeing : %ld\n",
             double(lower_bound) / scale_, disagreeing_global_indices_.size());

      max_lower_bound_ = std::max<double>(max_lower_bound_,double(lower_bound)/scale_);

      lower_bound_group_stats.addValue(lower_bound);
      if (lower_bound_group_stats.areGroupsPopulated()) {
        auto [first_group_max, second_group_max] =
            lower_bound_group_stats.getMaximums();
        if (second_group_max <= first_group_max) {
          std::cout << "breaking because max of this group's interval is less "
                       "than or qual to max in last last group's interval\n";
          opt_status = NO_FURTHER_PROGRESS;
          break;
        }
      }

      if (disagreeing_global_indices_.size() == 0) { // optimality condition
        std::cout << "breaking because of no disagreement\n";
        opt_status = OPTIMAL;
        break;
      }

      dual_cycle_list.addNode(getDualSolutionHash(
          disagreeing_global_indices_, lower_bound));
      if (dual_cycle_list.getMaxCycleCount() >
          max_cycle_count) { // at least one set of configurations likely
                             // repeated more than a specified number of times
        std::cout << "breaking because cycle detected\n";
        opt_status = NO_FURTHER_PROGRESS;
        break;
      }
    }
    printf(" === MAX === lower_bound : %ld\n", max_lower_bound);
    return opt_status;
  }

  template <int scale> void scaleProblem() {
    scale_ *= scale;
    for (auto &solver : solvers_) {
      thread_pool_.push([&] {
      solver->scaleProblem<scale>();
      });
    }
    thread_pool_.wait();
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      for (auto &constraint : constraints) {
        constraint.alpha *= scale;
      }
    }
  }

private:
  long getDualSolutionHash(const std::list<int> &disagreeing_global_indices,
                           long lower_bound) const {
    std::hash<long> hasher{};
    long hash = hasher(lower_bound);
    for (const auto &global_index : disagreeing_global_indices) {
      hash ^= hasher(global_index);
    }
    return hash;
  }

  std::list<int> runLagrangeMultipliersUpdateStep(int step_size,
                                                  bool use_momentum) {
    std::list<int> disagreeing_global_indices;
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      bool disagreement_exists = false;
      for (auto &constraint : constraints) {
        constraint.last_alpha =
            constraint.alpha; // record alpha before update
        int diff =
            solvers_[constraint.partition_index_target]->getMinCutSolution(
                constraint.local_index_target) -
            solvers_[constraint.partition_index_source]->getMinCutSolution(
                constraint.local_index_source);
        if (diff != 0) {
          if (use_momentum) {
            const double beta = .85;
            const int momentum_scale = 10;
            constraint.alpha_momentum =
                beta * constraint.alpha_momentum * beta + (1 - beta) * diff;
            constraint.alpha +=
                step_size *
                static_cast<int>(momentum_scale * constraint.alpha_momentum);
          } else {
            constraint.alpha += step_size * diff;
          }
          disagreement_exists = true;
        }
      }
      if (disagreement_exists) {
        disagreeing_global_indices.emplace_back(global_index);
      }
    }
    return std::move(disagreeing_global_indices);
  }

  void initializeDecomposition() {
    std::unordered_map</*global_index=*/int,
                       /*exists_in_partitions=*/std::set<int>>
        constrained_nodes;
    /**
     * step 0: parition graph into npartition_ partitions
     */
    std::vector<int> partitions_;
#ifdef HAVE_METIS
    partitions_ = metis_partition(npartition_,narc_,nnode_,arcs_);
#else 
    partitions_ = basic_graph_partition(npartition_,narc_,nnode_,arcs_);
#endif
    /**
     * step 1: distribute all arcs into one and only one sub graph
     */
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      if (s > t) {
        std::swap(s, t);
        std::swap(forward_capacity, backward_capacity);
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
      solvers_.emplace_back(std::make_unique<PrimalDualMinCutSolver>(std::move(min_cut_sub_graph.graph)));
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
      constraint_arc_map_.push_back({global_index,{}});
      auto &constraint_arcs = constraint_arc_map_.back().second;
      assert(partitions.size() >
             1); // requirement to be a proper constrained node
      for (const auto &partition_source :
           partitions) { // iterates in sorted order due to std::set
        for (const auto &partition_target :
             partitions) { // iterates in sorted order due to std::set
          if (partition_source >= partition_target) {
            continue;
          }
          int local_index_source =
              min_cut_sub_graphs_[partition_source].getNode(global_index);
          int local_index_target =
              min_cut_sub_graphs_[partition_target].getNode(global_index);
          constraint_arcs.emplace_back(
              /*alpha=*/0,
              /*last_alpha=*/0,
              /*alpha_momentum=*/0,
              /*partition_index_source=*/partition_source,
              /*partition_index_target=*/partition_target,
              /*local_index_source=*/local_index_source,
              /*local_index_target=*/local_index_target);
          auto arc_reference = --constraint_arcs.end();
          solvers_[partition_source]->addSourceDualDecompositionConstraint(
              arc_reference);
          solvers_[partition_target]->addTargetDualDecompositionConstraint(
              arc_reference);
        }
      }
      constrained_nodes_partition_counts.insert(partitions.size());
    }
    constraint_arc_map_.shrink_to_fit();
    printf("partition counts: ");
    for (const auto &count : constrained_nodes_partition_counts) {
      printf("%d,", count);
    }
    printf("\n");
    /**
     * step 5: create a min cut problem from the original problem to evaluate
     * the primal objective value
     */
    //primal_solver_ = std::make_unique<PrimalDualMinCutSolver>(nnode_,narc_,
    //   std::move(arcs_),
    //   std::move(arc_capacities_),
    //   std::move(terminal_capacities_));
    arcs_.clear();
    arcs_.shrink_to_fit();
    arc_capacities_.clear();
    arc_capacities_.shrink_to_fit();
    terminal_capacities_.clear();
    terminal_capacities_.shrink_to_fit();
  }

  /**
   * data passed into decomposition
   */
  int npartition_;
  int nnode_;
  int narc_;
  std::vector<int> arcs_;
  std::vector<int> arc_capacities_;
  std::vector<int> terminal_capacities_;

  /**
   * data structures needed for solving dual decomposition
   */
  std::vector<std::pair</*global_index=*/int,
                     std::list<DualDecompositionConstraintArc>>>
      constraint_arc_map_;

  std::vector<std::unique_ptr<PrimalDualMinCutSolver>> solvers_;
  std::unique_ptr<PrimalDualMinCutSolver>
      primal_solver_; // only used to evaluate primal value

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
  std::vector<bool> primal_solution_;
  long scale_;
  ThreadPool<void> thread_pool_;
  long solve_loop_time_;
  double max_lower_bound_;
  std::list<int> disagreeing_global_indices_;

  template <typename T> class ScalarStatisticsTracker {
  public:
    ScalarStatisticsTracker(size_t n) : n_(n), running_sum_(0), id_(0) {}

    void addValue(const T &value) {
      values_.push_back({value, id_});
      ordered_values_.insert(values_.back());
      id_ = (id_ + 1) % n_;
      running_sum_ += value;
      if (values_.size() > n_) {
        running_sum_ -= values_.front().first;
        ordered_values_.erase(values_.front());
        values_.pop_front();
      }
    }

    double getAverage() const { return static_cast<double>(running_sum_) / n_; }

    T getMaximum() const {
      if (!ordered_values_.size()) {
        return {};
      }
      return ordered_values_.rbegin()->first;
    }

  private:
    size_t n_;
    size_t id_;
    std::list<std::pair<T, size_t>> values_;
    T running_sum_;
    std::set<std::pair<T, size_t>> ordered_values_;
  };

  template <typename T> class TwoGroupScalarStatisticsTracker {
  public:
    TwoGroupScalarStatisticsTracker(size_t n)
        : n_(n), group_1_(n), group_2_(n), first_group(&group_1_),
          second_group(&group_2_), internal_counter_(0), is_ready_(false) {}

    void addValue(const T &value) {
      switch (internal_counter_ / n_) {
      case 0:
        first_group->addValue(value);
        break;
      case 1:
        second_group->addValue(value);
        break;
      default:
        assert(false);
      }
      internal_counter_++;
      if (internal_counter_ == 2 * n_) {
        std::swap(first_group, second_group);
        internal_counter_ = 0;
        is_ready_ = true;
      }
    }

    std::pair<T, T> getMaximums() const {
      return {second_group->getMaximum(), first_group->getMaximum()};
    }

    bool areGroupsPopulated() const {
      return is_ready_ && internal_counter_ == 0;
    }

  private:
    size_t n_;
    ScalarStatisticsTracker<T> group_1_, group_2_;
    ScalarStatisticsTracker<T> *first_group, *second_group;
    size_t internal_counter_;
    bool is_ready_;
  };
};

} // namespace mcpd3
