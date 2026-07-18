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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <list>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <numeric>

#include <measure/timer.h>
#include <decomp/constraint.h>
#include <decomp/halo_partition.h>
#include <decomp/lower_bound_certificate.h>
#include <decomp/optimization_schedule.h>
#include <decomp/partition_worker.h>
#include <graph/cycle.h>
#include <graph/partition.h>
#include <multithread/threadpool.h>
#include <primaldual/mcpd3.h>


namespace mcpd3 {

inline bool dualdecomp_progress_enabled() {
  const char *value = std::getenv("MCPD3_PROGRESS");
  return value != nullptr && value[0] != '\0' && std::string(value) != "0";
}

inline void dualdecomp_progress_report(
    const char *stage, long done, long total,
    std::chrono::steady_clock::time_point start) {
  if (!dualdecomp_progress_enabled()) {
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

inline void dualdecomp_progress_message(const std::string &message) {
  if (!dualdecomp_progress_enabled()) {
    return;
  }
  std::fprintf(stderr, "%s\n", message.c_str());
  std::fflush(stderr);
}

enum class DualDecompositionRegularizationScheme {
  SCALED_EPSILON,
  DISAGREEMENT_PLATEAU_EPSILON,
  NONE
};

struct DualDecompositionIterationRecord {
  long total_iteration = 0;
  int scale_iteration = 0;
  long objective_scale = 1;
  long step_size = 0;
  long effective_step_size = 0;
  Objective certified_lower_bound_raw = 0;
  Objective best_certified_lower_bound_raw = 0;
  Objective regularized_objective_raw = 0;
  long disagreement_count = 0;
  double disagreement_norm_sq = 0;
  Capacity regularization_strength = 0;
  Objective regularization_budget = 0;
  Objective regularization_contribution = 0;
  long solve_loop_microseconds = 0;
  long lagrange_update_microseconds = 0;
};

class DisagreementPlateauRegularizationTracker {
public:
  explicit DisagreementPlateauRegularizationTracker(int patience)
      : patience_(patience) {
    if (patience_ <= 0) {
      throw std::runtime_error("disagreement patience must be positive");
    }
    reset();
  }

  void reset() {
    active_ = false;
    has_observation_ = false;
    has_regularization_pulse_ = false;
    best_disagreement_count_ = std::numeric_limits<long>::max();
    last_improvement_iteration_ = 0;
    last_regularization_pulse_iteration_ = 0;
  }

  bool observe(int iteration, long disagreement_count) {
    if (iteration < 0 || disagreement_count < 0) {
      throw std::runtime_error(
          "disagreement plateau observations must be non-negative");
    }
    if (!has_observation_ ||
        disagreement_count < best_disagreement_count_) {
      has_observation_ = true;
      best_disagreement_count_ = disagreement_count;
      last_improvement_iteration_ = iteration;
      return false;
    }
    const int window_start =
        has_regularization_pulse_
            ? std::max(last_improvement_iteration_,
                       last_regularization_pulse_iteration_)
            : last_improvement_iteration_;
    if (iteration - window_start >= patience_) {
      active_ = true;
      has_regularization_pulse_ = true;
      last_regularization_pulse_iteration_ = iteration;
      return true;
    }
    return false;
  }

  bool active() const { return active_; }
  long bestDisagreementCount() const { return best_disagreement_count_; }
  int iterationsSinceImprovement(int iteration) const {
    return has_observation_ ? iteration - last_improvement_iteration_ : 0;
  }

private:
  int patience_;
  bool active_ = false;
  bool has_observation_ = false;
  bool has_regularization_pulse_ = false;
  long best_disagreement_count_ = std::numeric_limits<long>::max();
  int last_improvement_iteration_ = 0;
  int last_regularization_pulse_iteration_ = 0;
};

struct DualDecompositionOptions {
  int num_optimization_scales = 5;
  int max_iteration_count = 10000;
  long max_total_iteration_count = 0;
  int max_cycle_count = 2;
  long initial_step_size = 10000;
  int patience = 10;
  int disagreement_patience = 10;
  bool legacy_patience = false;
  bool exhaust_scale_iterations = false;
  bool exhaust_regularized_scale_iterations = false;
  bool use_momentum = true;
  bool enable_group_stopping = true;
  bool track_primal_upper_bound = true;
  bool emit_partition_packages = true;
  bool construct_solvers = true;
  bool materialize_all_partition_nodes = false;
  bool saturate_capacity_overflow = false;
  bool verbose = true;
  long min_step_size = 1;
  long objective_scale = 1;
  size_t thread_count = 0;
  DualDecompositionRegularizationScheme regularization_scheme =
      DualDecompositionRegularizationScheme::SCALED_EPSILON;
  int scaled_epsilon_max_step_size = 10;
  int scaled_epsilon_strength_cap = 0;
  Objective regularization_budget_limit = 0;
  bool promote_objective_scale_on_overbudget = true;
  int max_objective_scale_promotions = 4;
  bool retry_unit_step_without_momentum = false;
  bool randomize_initial_alphas = false;
  long initial_alpha_random_radius = 0;
  unsigned int initial_alpha_random_seed = 0;
  CanonicalCutSelection canonical_cut_selection =
      CanonicalCutSelection::SOLVER_DEFAULT;
  bool force_full_mincut_recompute = false;
  bool track_arc_flow_updates = false;
  int halo_depth = 1;
  std::vector<std::uint64_t> partition_edge_weights;
  std::vector<int> partition_labels;
  std::vector<int> reference_cut_labels;
  ReferenceCutSelection reference_cut_selection =
      ReferenceCutSelection::CLOSEST_EXACT;
  long reference_cut_check_interval = 1;
  std::function<void(const DualDecompositionIterationRecord &)>
      iteration_callback;
};

class DualDecomposition {
public:
  DualDecomposition(int npartition, int nnode, int narc, std::vector<int> arcs,
                    std::vector<Capacity> arc_capacities,
                    std::vector<Capacity> terminal_capacities,
                    DualDecompositionOptions options = {})
      : npartition_(npartition), nnode_(nnode), narc_(narc),
        arcs_(std::move(arcs)), arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)),
        original_arcs_(options.track_primal_upper_bound ? arcs_
                                                        : std::vector<int>()),
        original_arc_capacities_(options.track_primal_upper_bound
                                     ? arc_capacities_
                                     : std::vector<Capacity>()),
        original_terminal_capacities_(options.track_primal_upper_bound
                                          ? terminal_capacities_
                                          : std::vector<Capacity>()),
        min_cut_sub_graphs_(npartition_),
        partition_packages_(options.emit_partition_packages ? npartition_ : 0),
        primal_solution_(nnode_), scale_(1),
        options_(options),
        thread_pool_(
            resolveThreadCount(npartition_, options_.thread_count)),
        solve_loop_time_(0),
        lagrange_update_time_(0),
        max_lower_bound_(std::numeric_limits<double>::lowest()),
        max_lower_bound_raw_(0), max_regularized_objective_raw_(0),
        best_upper_bound_(0), current_upper_bound_(0),
        has_max_lower_bound_raw_(false),
        has_max_regularized_objective_raw_(false),
        has_best_upper_bound_(false), has_current_upper_bound_(false),
        last_original_objective_raw_(0),
        last_certified_lower_bound_raw_(0),
        last_regularized_objective_raw_(0),
        last_disagreement_count_(0),
        last_disagreement_norm_sq_(0),
        last_regularization_budget_(0),
        last_regularization_contribution_(0),
        last_regularization_anchor_sink_count_(0),
        last_regularization_active_sink_count_(0),
        total_optimization_iterations_(0),
        objective_scale_promotion_count_(0),
        halo_objective_multiplier_(1),
        warned_regularization_budget_exceeded_(false) {
    validateOptions();
    initializeDecomposition();
  }

  template <typename InputCapacity,
            std::enable_if_t<!std::is_same_v<InputCapacity, Capacity>, int> = 0>
  DualDecomposition(int npartition, int nnode, int narc,
                    std::vector<int> arcs,
                    const std::vector<InputCapacity> &arc_capacities,
                    const std::vector<InputCapacity> &terminal_capacities,
                    DualDecompositionOptions options = {})
      : DualDecomposition(npartition, nnode, narc, std::move(arcs),
                          capacity_vector_from(arc_capacities),
                          capacity_vector_from(terminal_capacities), options) {}

  DualDecomposition(int npartition, MinCutGraph min_cut_graph,
                    DualDecompositionOptions options = {})
      : DualDecomposition(npartition, min_cut_graph.nnode, min_cut_graph.narc,
                          std::move(min_cut_graph.arcs),
                          std::move(min_cut_graph.arc_capacities),
                          std::move(min_cut_graph.terminal_capacities),
                          options) {}

  long getTotalSolveLoopTime() const { return solve_loop_time_; }
  long getTotalLagrangeUpdateTime() const { return lagrange_update_time_; }
  long getScale() const { return scale_; }
  double getBestLowerBound() const { return max_lower_bound_; }
  Objective getBestLowerBoundRaw() const { return max_lower_bound_raw_; }
  double getBestCertifiedLowerBound() const { return max_lower_bound_; }
  Objective getBestCertifiedLowerBoundRaw() const {
    return max_lower_bound_raw_;
  }
  double getBestRegularizedObjective() const {
    return !has_max_regularized_objective_raw_
               ? -std::numeric_limits<double>::infinity()
               : integer_to_double(max_regularized_objective_raw_) / scale_;
  }
  Objective getBestRegularizedObjectiveRaw() const {
    return max_regularized_objective_raw_;
  }
  Objective getBestUpperBoundRaw() const { return best_upper_bound_; }
  bool hasBestUpperBound() const { return has_best_upper_bound_; }
  double getBestUpperBound() const {
    return !has_best_upper_bound_
               ? std::numeric_limits<double>::infinity()
               : integer_to_double(best_upper_bound_) / scale_;
  }
  Objective getCurrentUpperBoundRaw() const { return current_upper_bound_; }
  bool hasCurrentUpperBound() const { return has_current_upper_bound_; }
  Objective getLastOriginalObjectiveRaw() const {
    return last_original_objective_raw_;
  }
  Objective getLastCertifiedLowerBoundRaw() const {
    return last_certified_lower_bound_raw_;
  }
  Objective getLastRegularizedObjectiveRaw() const {
    return last_regularized_objective_raw_;
  }
  long getLastDisagreementCount() const { return last_disagreement_count_; }
  double getLastDisagreementNormSq() const {
    return last_disagreement_norm_sq_;
  }
  Objective getLastRegularizationBudget() const {
    return last_regularization_budget_;
  }
  Objective getLastRegularizationContribution() const {
    return last_regularization_contribution_;
  }
  long getLastRegularizationAnchorSinkCount() const {
    return last_regularization_anchor_sink_count_;
  }
  long getLastRegularizationActiveSinkCount() const {
    return last_regularization_active_sink_count_;
  }
  long getTotalOptimizationIterations() const {
    return total_optimization_iterations_;
  }
  long getReferenceDecodeCount() const {
    long count = 0;
    for (const auto &solver : solvers_) {
      count += solver->getReferenceDecodeCount();
    }
    return count;
  }
  long getReferenceCurrentCutHitCount() const {
    long count = 0;
    for (const auto &solver : solvers_) {
      count += solver->getReferenceCurrentCutHitCount();
    }
    return count;
  }
  long getReferenceExactHitCount() const {
    long count = 0;
    for (const auto &solver : solvers_) {
      count += solver->getReferenceExactHitCount();
    }
    return count;
  }
  long getReferenceClosureCount() const {
    long count = 0;
    for (const auto &solver : solvers_) {
      count += solver->getReferenceClosureCount();
    }
    return count;
  }
  long getReferenceDecodeTimeMicroseconds() const {
    long elapsed_us = 0;
    for (const auto &solver : solvers_) {
      elapsed_us += solver->getReferenceDecodeTimeMicroseconds();
    }
    return elapsed_us;
  }
  long getObjectiveScalePromotionCount() const {
    return objective_scale_promotion_count_;
  }
  long getUnitStepNoMomentumRetryCount() const {
    return unit_step_no_momentum_retry_count_;
  }
  long getHaloObjectiveMultiplier() const {
    return halo_objective_multiplier_;
  }
  const std::vector<int> &getPartitionLabels() const {
    return partition_labels_.empty() ? options_.partition_labels
                                     : partition_labels_;
  }
  std::vector<std::uint64_t> getArcFlowUpdateCounts() const {
    requireConstructedSolvers("getArcFlowUpdateCounts");
    std::vector<std::uint64_t> counts(static_cast<size_t>(narc_), 0);
    for (int arc = 0; arc < narc_; ++arc) {
      if (options_.halo_depth == 1) {
        const auto &location = arc_locations_[static_cast<size_t>(arc)];
        const auto &local_counts =
            solvers_[static_cast<size_t>(location.partition)]
                ->getArcFlowUpdateCounts();
        counts[static_cast<size_t>(arc)] =
            local_counts[static_cast<size_t>(location.local_arc)];
        continue;
      }
      for (const auto &location :
           halo_arc_locations_[static_cast<size_t>(arc)]) {
        const auto &local_counts =
            solvers_[static_cast<size_t>(location.partition)]
                ->getArcFlowUpdateCounts();
        const std::uint64_t local_count =
            local_counts[static_cast<size_t>(location.local_arc)];
        auto &total = counts[static_cast<size_t>(arc)];
        if (local_count > std::numeric_limits<std::uint64_t>::max() - total) {
          throw std::overflow_error("halo arc flow-update count overflow");
        }
        total += local_count;
      }
    }
    return counts;
  }
  void resetArcFlowUpdateCounts() {
    requireConstructedSolvers("resetArcFlowUpdateCounts");
    for (auto &solver : solvers_) {
      solver->resetArcFlowUpdateCounts();
    }
  }
  int getConfiguredNumOptimizationScales() const {
    return options_.num_optimization_scales;
  }
  long getConfiguredInitialStepSize() const {
    return options_.initial_step_size;
  }
  bool getConfiguredExhaustRegularizedScaleIterations() const {
    return options_.exhaust_regularized_scale_iterations;
  }
  void configureOptimizationSchedule(
      int num_optimization_scales, long initial_step_size,
      bool exhaust_regularized_scale_iterations) {
    if (num_optimization_scales <= 0) {
      throw std::runtime_error(
          "optimization scale count must be positive");
    }
    if (initial_step_size <= 0) {
      throw std::runtime_error("initial step size must be positive");
    }
    options_.num_optimization_scales = num_optimization_scales;
    options_.initial_step_size = initial_step_size;
    options_.exhaust_regularized_scale_iterations =
        exhaust_regularized_scale_iterations;
  }
  const std::vector<PartitionPackage> &getPartitionPackages() const {
    if (!options_.emit_partition_packages) {
      throw std::runtime_error(
          "partition package export is disabled for this DualDecomposition");
    }
    return partition_packages_;
  }

  std::vector<DualDecompositionConstraintSnapshot>
  getConstraintSnapshots() const {
    std::vector<DualDecompositionConstraintSnapshot> snapshots;
    int constraint_id = 0;
    for (const auto &[global_index, constraints] : constraint_arc_map_) {
      for (const auto &constraint : constraints) {
        snapshots.push_back(DualDecompositionConstraintSnapshot{
            /*constraint_id=*/constraint_id++,
            /*global_node_id=*/global_index,
            /*partition_index_source=*/constraint.partition_index_source,
            /*partition_index_target=*/constraint.partition_index_target,
            /*local_index_source=*/constraint.local_index_source,
            /*local_index_target=*/constraint.local_index_target,
            /*alpha=*/constraint.alpha,
            /*last_alpha=*/constraint.last_alpha,
            /*alpha_momentum=*/constraint.alpha_momentum});
      }
    }
    return snapshots;
  }

  std::vector<DualDecompositionPartitionSnapshot>
  getPartitionSnapshots() const {
    requireConstructedSolvers("getPartitionSnapshots");
    std::vector<DualDecompositionPartitionSnapshot> snapshots;
    snapshots.reserve(solvers_.size());
    for (size_t partition_id = 0; partition_id < solvers_.size();
         ++partition_id) {
      const auto &solver = solvers_[partition_id];
      const int local_node_count =
          min_cut_sub_graphs_[partition_id].graph.nnode;
      DualDecompositionPartitionSnapshot snapshot;
      snapshot.partition_id = static_cast<int>(partition_id);
      snapshot.lower_bound = solver->getMinCutValue();
      snapshot.regularization_budget =
          solver->getLastRegularizationBudget();
      snapshot.regularization_contribution =
          solver->getLastRegularizationContribution();
      snapshot.regularization_anchor_sink_count =
          solver->getLastRegularizationAnchorSinkCount();
      snapshot.regularization_active_sink_count =
          solver->getLastRegularizationActiveSinkCount();
      snapshot.local_labels.reserve(static_cast<size_t>(local_node_count));
      for (int local_index = 0; local_index < local_node_count; ++local_index) {
        snapshot.local_labels.push_back(
            solver->getMinCutSolution(local_index));
      }
      snapshots.push_back(std::move(snapshot));
    }
    return snapshots;
  }

  struct FlowWarmStart {
    std::vector<PrimalDualMinCutSolver::FlowWarmStart> partitions;
    std::vector<DualDecompositionConstraintSnapshot> constraints;
  };

  FlowWarmStart captureFlowWarmStart() const {
    requireConstructedSolvers("captureFlowWarmStart");
    FlowWarmStart state;
    state.partitions.reserve(solvers_.size());
    for (const auto &solver : solvers_) {
      state.partitions.push_back(solver->captureFlowWarmStart());
    }
    state.constraints = getConstraintSnapshots();
    return state;
  }

  void restoreFlowWarmStart(const FlowWarmStart &state) {
    requireConstructedSolvers("restoreFlowWarmStart");
    if (state.partitions.size() != solvers_.size()) {
      throw std::runtime_error("DD flow warm start partition count mismatch");
    }

    size_t constraint_index = 0;
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      for (auto &constraint : constraints) {
        if (constraint_index >= state.constraints.size()) {
          throw std::runtime_error(
              "DD flow warm start constraint count mismatch");
        }
        const auto &snapshot = state.constraints[constraint_index];
        if (snapshot.constraint_id != static_cast<int>(constraint_index) ||
            snapshot.global_node_id != global_index ||
            snapshot.partition_index_source !=
                constraint.partition_index_source ||
            snapshot.partition_index_target !=
                constraint.partition_index_target ||
            snapshot.local_index_source != constraint.local_index_source ||
            snapshot.local_index_target != constraint.local_index_target) {
          throw std::runtime_error(
              "DD flow warm start constraint topology mismatch");
        }
        constraint.alpha = snapshot.alpha;
        constraint.last_alpha = snapshot.last_alpha;
        constraint.alpha_momentum = snapshot.alpha_momentum;
        ++constraint_index;
      }
    }
    if (constraint_index != state.constraints.size()) {
      throw std::runtime_error("DD flow warm start constraint count mismatch");
    }
    for (size_t partition = 0; partition < solvers_.size(); ++partition) {
      solvers_[partition]->restoreFlowWarmStart(state.partitions[partition]);
    }
  }

  void replaceProblemCapacities(
      const std::vector<Capacity> &arc_capacities,
      const std::vector<Capacity> &terminal_capacities,
      bool preserve_alpha_state = true,
      bool preserve_flow_state = true,
      const Objective &flow_scale_numerator = 1,
      const Objective &flow_scale_denominator = 1) {
    requireConstructedSolvers("replaceProblemCapacities");
    if (arc_capacities.size() != static_cast<size_t>(2 * narc_)) {
      throw std::runtime_error("replacement arc capacity count mismatch");
    }
    if (terminal_capacities.size() != static_cast<size_t>(nnode_)) {
      throw std::runtime_error(
          "replacement terminal capacity count mismatch");
    }
    for (const Capacity &capacity : arc_capacities) {
      if (capacity < 0) {
        throw std::runtime_error(
            "replacement arc capacities must be non-negative");
      }
    }

    std::vector<std::vector<Capacity>> local_arc_capacities(solvers_.size());
    std::vector<std::vector<Capacity>> local_terminal_capacities(
        solvers_.size());
    for (size_t partition = 0; partition < solvers_.size(); ++partition) {
      local_arc_capacities[partition].assign(
          static_cast<size_t>(2 * local_arc_counts_[partition]), 0);
      local_terminal_capacities[partition].assign(
          static_cast<size_t>(local_node_counts_[partition]), 0);
    }
    for (int arc = 0; arc < narc_; ++arc) {
      const Capacity input_forward = arc_capacities[2 * arc];
      const Capacity input_backward = arc_capacities[2 * arc + 1];
      if (options_.halo_depth == 1) {
        const ArcLocation &location = arc_locations_[arc];
        auto &local = local_arc_capacities[location.partition];
        local[2 * location.local_arc] =
            location.swapped ? input_backward : input_forward;
        local[2 * location.local_arc + 1] =
            location.swapped ? input_forward : input_backward;
        continue;
      }
      const long objective_factor =
          halo_arc_objective_factors_[static_cast<size_t>(arc)];
      for (const ArcLocation &location : halo_arc_locations_[arc]) {
        const Capacity local_forward = checked_scale_capacity(
            input_forward, objective_factor,
            options_.saturate_capacity_overflow);
        const Capacity local_backward = checked_scale_capacity(
            input_backward, objective_factor,
            options_.saturate_capacity_overflow);
        auto &local = local_arc_capacities[location.partition];
        local[2 * location.local_arc] =
            location.swapped ? local_backward : local_forward;
        local[2 * location.local_arc + 1] =
            location.swapped ? local_forward : local_backward;
      }
    }
    for (int node = 0; node < nnode_; ++node) {
      if (options_.halo_depth == 1) {
        const TerminalLocation &location = terminal_locations_[node];
        if (location.partition < 0 || location.local_node < 0) {
          if (terminal_capacities[node] != 0) {
            throw std::runtime_error(
                "replacement terminal activates an absent isolated node");
          }
          continue;
        }
        local_terminal_capacities[location.partition][location.local_node] =
            terminal_capacities[node];
        continue;
      }
      const auto &locations = halo_terminal_locations_[node];
      if (locations.empty()) {
        if (terminal_capacities[node] != 0) {
          throw std::runtime_error(
              "replacement terminal activates an absent isolated node");
        }
        continue;
      }
      const long objective_factor =
          halo_terminal_objective_factors_[static_cast<size_t>(node)];
      for (const TerminalLocation &location : locations) {
        local_terminal_capacities[location.partition][location.local_node] =
            checked_scale_capacity(terminal_capacities[node], objective_factor,
                                   options_.saturate_capacity_overflow);
      }
    }

    if (preserve_flow_state &&
        flow_scale_numerator != flow_scale_denominator) {
      for (size_t partition = 0; partition < solvers_.size(); ++partition) {
        solvers_[partition]->validateProblemCapacityReplacement(
            local_arc_capacities[partition],
            local_terminal_capacities[partition], preserve_flow_state,
            flow_scale_numerator, flow_scale_denominator);
      }
    }

    for (size_t partition = 0; partition < solvers_.size(); ++partition) {
      solvers_[partition]->replaceProblemCapacities(
          local_arc_capacities[partition],
          local_terminal_capacities[partition], preserve_flow_state,
          flow_scale_numerator, flow_scale_denominator);
      if (options_.emit_partition_packages) {
        partition_packages_[partition].arc_capacities =
            local_arc_capacities[partition];
        partition_packages_[partition].terminal_capacities =
            local_terminal_capacities[partition];
      }
    }

    if (!preserve_alpha_state) {
      for (auto &[global_index, constraints] : constraint_arc_map_) {
        (void)global_index;
        for (auto &constraint : constraints) {
          constraint.alpha = 0;
          constraint.last_alpha = 0;
          constraint.alpha_momentum = 0;
        }
      }
    }

    if (options_.track_primal_upper_bound) {
      original_arc_capacities_.clear();
      original_arc_capacities_.reserve(arc_capacities.size());
      for (const Capacity &capacity : arc_capacities) {
        original_arc_capacities_.push_back(checked_scale_capacity(
            capacity, halo_objective_multiplier_,
            options_.saturate_capacity_overflow));
      }
      original_terminal_capacities_.clear();
      original_terminal_capacities_.reserve(terminal_capacities.size());
      for (const Capacity &capacity : terminal_capacities) {
        original_terminal_capacities_.push_back(checked_scale_capacity(
            capacity, halo_objective_multiplier_,
            options_.saturate_capacity_overflow));
      }
    }
    solve_loop_time_ = 0;
    lagrange_update_time_ = 0;
    max_lower_bound_ = std::numeric_limits<double>::lowest();
    max_lower_bound_raw_ = 0;
    max_regularized_objective_raw_ = 0;
    best_upper_bound_ = 0;
    current_upper_bound_ = 0;
    has_max_lower_bound_raw_ = false;
    has_max_regularized_objective_raw_ = false;
    has_best_upper_bound_ = false;
    has_current_upper_bound_ = false;
    last_original_objective_raw_ = 0;
    last_certified_lower_bound_raw_ = 0;
    last_regularized_objective_raw_ = 0;
    last_disagreement_count_ = 0;
    last_disagreement_norm_sq_ = 0;
    last_regularization_budget_ = 0;
    last_regularization_contribution_ = 0;
    last_regularization_anchor_sink_count_ = 0;
    last_regularization_active_sink_count_ = 0;
    total_optimization_iterations_ = 0;
    objective_scale_promotion_count_ = 0;
    warned_regularization_budget_exceeded_ = false;
    disagreeing_global_indices_.clear();
  }

  template <typename InputCapacity,
            std::enable_if_t<!std::is_same_v<InputCapacity, Capacity>, int> = 0>
  void replaceProblemCapacities(
      const std::vector<InputCapacity> &arc_capacities,
      const std::vector<InputCapacity> &terminal_capacities,
      bool preserve_alpha_state = true, bool preserve_flow_state = true,
      const Objective &flow_scale_numerator = 1,
      const Objective &flow_scale_denominator = 1) {
    replaceProblemCapacities(capacity_vector_from(arc_capacities),
                             capacity_vector_from(terminal_capacities),
                             preserve_alpha_state, preserve_flow_state,
                             flow_scale_numerator, flow_scale_denominator);
  }

  int regularizationStrengthForStepSize(long step_size) const {
    if (options_.regularization_scheme !=
        DualDecompositionRegularizationScheme::SCALED_EPSILON) {
      return 0;
    }
    if (step_size > options_.scaled_epsilon_max_step_size) {
      return 0;
    }
    const int step_strength = static_cast<int>(step_size);
    return options_.scaled_epsilon_strength_cap > 0
               ? std::min(step_strength, options_.scaled_epsilon_strength_cap)
               : step_strength;
  }

  int plateauRegularizationStrength() const { return 1; }

  long getDisagreementPlateauActivationCount() const {
    return disagreement_plateau_activation_count_;
  }

  void runPrimalSolutionDecodingStep(bool do_narrow_band_decode = false) {
    requireConstructedSolvers("primal decoding");
    for (int i = 0; i < npartition_; ++i) {
      const auto &min_cut_sub_graph = min_cut_sub_graphs_[i];
      const auto &solver = solvers_[i];
      for (int local_index = 0;
           local_index < static_cast<int>(min_cut_sub_graph.local_to_global.size());
           ++local_index) {
        const int global_index = min_cut_sub_graph.local_to_global[local_index];
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
      std::cout << "recalc primal: "
                << integer_to_string(primal_solver_->getMinCutValue())
                << "\n";
    }
  }

  enum OptimizationStatus {
    OPTIMAL,
    NO_FURTHER_PROGRESS,
    ITERATION_COUNT_EXCEEDED,
    REGULARIZATION_BUDGET_EXCEEDED
  };

  template <bool attempt_decoding, typename Decoder>
  void solve(Decoder decoder) {
    requireConstructedSolvers("solve");
    long step_size = options_.initial_step_size;
    scale_ = options_.objective_scale;
    total_optimization_iterations_ = 0;
    unit_step_no_momentum_retry_count_ = 0;
    disagreement_plateau_activation_count_ = 0;
    int iscale = 0;
    bool unit_step_no_momentum_retry_attempted = false;
    int schedule_level_count = options_.num_optimization_scales;
    while (iscale < schedule_level_count && step_size >= 1) {
      OptimizationStatus status;
      auto run_opt_scale_time = time_lambda([&] {
        status = runOptimizationScale(options_.max_iteration_count, step_size,
                                      options_.max_cycle_count,
                                      options_.use_momentum);
      });
      if (options_.verbose) {
        printf("run optimization scale time: %lums\n",
               run_opt_scale_time.count());
      }
      if (totalIterationBudgetExhausted()) {
        break;
      }
      if (status == REGULARIZATION_BUDGET_EXCEEDED &&
          tryPromoteObjectiveScale(/*factor=*/10, &step_size)) {
        schedule_level_count = std::max(
            schedule_level_count, optimizationScheduleLevelCount(step_size));
        iscale = 0;
        unit_step_no_momentum_retry_attempted = false;
        continue;
      }
      if (status == mcpd3::DualDecomposition::OPTIMAL) {
        break;
      }
      if (attempt_decoding && iscale > 0) { // decoding requested
        runPrimalSolutionDecodingStep();
        if (decoder(primal_solution_, max_lower_bound_,
                    disagreeing_global_indices_)) {
          break;
        }
      }
      if (step_size == 1 && options_.use_momentum &&
          options_.retry_unit_step_without_momentum &&
          !unit_step_no_momentum_retry_attempted) {
        unit_step_no_momentum_retry_attempted = true;
        ++unit_step_no_momentum_retry_count_;
        if (dualdecomp_progress_enabled()) {
          std::fprintf(stderr,
                       "mcpd3_progress stage=dd_unit_step_no_momentum_retry "
                       "scale=%ld retry_count=%ld\n",
                       scale_, unit_step_no_momentum_retry_count_);
          std::fflush(stderr);
        }
        status = runOptimizationScale(options_.max_iteration_count,
                                      /*step_size=*/1,
                                      options_.max_cycle_count,
                                      /*use_momentum=*/false);
        if (status == mcpd3::DualDecomposition::OPTIMAL) {
          break;
        }
        if (status == REGULARIZATION_BUDGET_EXCEEDED &&
            tryPromoteObjectiveScale(/*factor=*/10, &step_size)) {
          schedule_level_count = std::max(
              schedule_level_count, optimizationScheduleLevelCount(step_size));
          iscale = 0;
          unit_step_no_momentum_retry_attempted = false;
          continue;
        }
      }
      if (step_size == 1 &&
          tryPromoteObjectiveScale(/*factor=*/10, &step_size)) {
        schedule_level_count = std::max(
            schedule_level_count, optimizationScheduleLevelCount(step_size));
        iscale = 0;
        unit_step_no_momentum_retry_attempted = false;
        continue;
      }
      step_size = nextOptimizationScheduleValue(step_size);
      ++iscale;
    }
  }

  void solve() {
    auto null_decoder =
        [=](const std::vector<bool> &cut, double max_lower_bound,
            const std::list<int> &disagreeing_global_indices) -> bool {
      return false;
    };
    solve<false>(null_decoder);
  }

  struct LagrangeUpdateStats {
    std::list<int> disagreeing_global_indices;
    long disagreement_count = 0;
    double disagreement_norm_sq = 0;
    long effective_step_size = 0;
  };

  OptimizationStatus runOptimizationScale(int nstep, long step_size,
                                          int max_cycle_count = 2,
                                          int use_momentum = false) {
    requireConstructedSolvers("runOptimizationScale");
    OptimizationStatus opt_status = ITERATION_COUNT_EXCEEDED;
    const bool report_progress = dualdecomp_progress_enabled();
    const auto scale_start = std::chrono::steady_clock::now();
    const int num_stats_in_group = 10;
    TwoGroupScalarStatisticsTracker<double> lower_bound_group_stats(
        num_stats_in_group);
    CycleCountingList dual_cycle_list;
    Objective max_lower_bound = 0;
    bool has_scale_max_lower_bound = false;
    int last_improvement_iter = 0;
    const bool disagreement_plateau_mode =
        options_.regularization_scheme ==
        DualDecompositionRegularizationScheme::
            DISAGREEMENT_PLATEAU_EPSILON;
    std::unique_ptr<DisagreementPlateauRegularizationTracker>
        disagreement_plateau_tracker;
    if (disagreement_plateau_mode) {
      disagreement_plateau_tracker =
          std::make_unique<DisagreementPlateauRegularizationTracker>(
              options_.disagreement_patience);
    }
    int regularization_strength =
        regularizationStrengthForStepSize(step_size);
    for (auto &solver_uptr : solvers_) {
      solver_uptr->setRegularizationStrength(regularization_strength);
    }
    for (int i = 0; i < nstep; ++i) {
      if (totalIterationBudgetExhausted()) {
        break;
      }
      ++total_optimization_iterations_;

      std::vector<Objective> lower_bound_terms(solvers_.size(), 0);
      std::vector<Objective> regularization_budget_terms(solvers_.size(), 0);
      std::vector<Objective> regularization_contribution_terms(solvers_.size(), 0);
      std::vector<long> regularization_anchor_count_terms(solvers_.size(), 0);
      std::vector<long> regularization_active_count_terms(solvers_.size(), 0);
      auto solve_loop_time = time_lambda([&] {
        for (size_t solver_index = 0; solver_index < solvers_.size();
             ++solver_index) {
          auto *solver = solvers_[solver_index].get();
          auto *lower_result = &lower_bound_terms[solver_index];
          auto *regularization_budget_result =
              &regularization_budget_terms[solver_index];
          auto *regularization_contribution_result =
              &regularization_contribution_terms[solver_index];
          auto *regularization_anchor_count_result =
              &regularization_anchor_count_terms[solver_index];
          auto *regularization_active_count_result =
              &regularization_active_count_terms[solver_index];
          thread_pool_.push([solver, lower_result,
                             regularization_budget_result,
                             regularization_contribution_result,
                             regularization_anchor_count_result,
                             regularization_active_count_result] {
            solver->solve();
            *lower_result = solver->getMinCutValue();
            *regularization_budget_result =
                solver->getLastRegularizationBudget();
            *regularization_contribution_result =
                solver->getLastRegularizationContribution();
            *regularization_anchor_count_result =
                solver->getLastRegularizationAnchorSinkCount();
            *regularization_active_count_result =
                solver->getLastRegularizationActiveSinkCount();
          });
        }
        thread_pool_.wait();
      });
      solve_loop_time_ += solve_loop_time.count();
      Objective original_objective = 0;
      last_regularization_budget_ = 0;
      last_regularization_contribution_ = 0;
      for (const auto &term : lower_bound_terms) {
        original_objective = checked_add(
            original_objective, term, "local objective sum overflow");
      }
      for (const auto &term : regularization_budget_terms) {
        last_regularization_budget_ = checked_add(
            last_regularization_budget_, term,
            "regularization budget sum overflow");
      }
      for (const auto &term : regularization_contribution_terms) {
        last_regularization_contribution_ = checked_add(
            last_regularization_contribution_, term,
            "regularization contribution sum overflow");
      }
      last_regularization_anchor_sink_count_ =
          std::accumulate(regularization_anchor_count_terms.begin(),
                          regularization_anchor_count_terms.end(),
                          static_cast<long>(0));
      last_regularization_active_sink_count_ =
          std::accumulate(regularization_active_count_terms.begin(),
                          regularization_active_count_terms.end(),
                          static_cast<long>(0));
      const Objective regularized_objective =
          regularizedObjectiveRaw(original_objective,
                                  last_regularization_contribution_);
      Objective lower_bound = certifiedOriginalLowerBoundRaw(
          original_objective, last_regularization_contribution_,
          last_regularization_budget_);
      last_original_objective_raw_ = original_objective;
      last_certified_lower_bound_raw_ = lower_bound;
      last_regularized_objective_raw_ = regularized_objective;
      warnIfRegularizationBudgetExceeded(last_regularization_budget_,
                                         regularization_strength);
      const bool regularization_budget_exceeded =
          isRegularizationBudgetExceeded(last_regularization_budget_,
                                         regularization_strength);
      if (regularization_budget_exceeded) {
        if (report_progress) {
          std::fprintf(
              stderr,
              "mcpd3_progress stage=dd_solve_stop "
              "reason=regularization_budget_exceeded iter=%d "
              "budget=%s limit=%s step_size=%ld scale=%ld\n",
              i, integer_to_string(last_regularization_budget_).c_str(),
              integer_to_string(regularizationBudgetLimit()).c_str(),
              step_size, scale_);
          std::fflush(stderr);
        }
        opt_status = REGULARIZATION_BUDGET_EXCEEDED;
        break;
      }
      if (options_.track_primal_upper_bound) {
        current_upper_bound_ = updatePrimalUpperBound();
        has_current_upper_bound_ = true;
      }

      LagrangeUpdateStats update_stats;
      auto lagrange_update_time = time_lambda([&] {
        update_stats =
            runLagrangeMultipliersUpdateStep(step_size, use_momentum,
                                             lower_bound);
        disagreeing_global_indices_ =
            std::move(update_stats.disagreeing_global_indices);
          });
      lagrange_update_time_ += lagrange_update_time.count();
      last_disagreement_count_ = update_stats.disagreement_count;
      last_disagreement_norm_sq_ = update_stats.disagreement_norm_sq;
      if (update_stats.disagreement_count == 0 &&
          last_regularization_budget_ < Objective(options_.objective_scale)) {
        lower_bound = original_objective;
        last_certified_lower_bound_raw_ = lower_bound;
      }

      const bool plateau_was_active =
          disagreement_plateau_mode && disagreement_plateau_tracker->active();
      const bool plateau_regularization_pulse_now =
          disagreement_plateau_mode && update_stats.disagreement_count > 0 &&
          disagreement_plateau_tracker->observe(
              i, update_stats.disagreement_count);
      const bool plateau_activated_now =
          plateau_regularization_pulse_now && !plateau_was_active;
      if (plateau_regularization_pulse_now) {
        regularization_strength = plateauRegularizationStrength();
        for (auto &solver_uptr : solvers_) {
          solver_uptr->setRegularizationStrength(regularization_strength);
        }
        last_improvement_iter = i;
        if (plateau_activated_now) {
          ++disagreement_plateau_activation_count_;
        }
        if (report_progress) {
          std::fprintf(
              stderr,
              "mcpd3_progress stage=%s "
              "scale=%ld step_size=%ld iter=%d disagreement_count=%ld "
              "disagreement_patience=%d regularization_strength=%s\n",
              plateau_activated_now ? "dd_disagreement_plateau_activate"
                                    : "dd_disagreement_plateau_pulse",
              scale_, step_size, i, update_stats.disagreement_count,
              options_.disagreement_patience,
              integer_to_string(
                  widen_capacity(regularization_strength)).c_str());
          std::fflush(stderr);
        }
      }
      const bool standard_early_exit_enabled =
          !disagreement_plateau_mode ||
          disagreement_plateau_tracker->active();
      const Objective best_lower_bound =
          !has_scale_max_lower_bound || lower_bound > max_lower_bound
              ? lower_bound
              : max_lower_bound;
      if (report_progress) {
        const Objective best_regularized_objective =
            !has_max_regularized_objective_raw_ ||
                    regularized_objective > max_regularized_objective_raw_
                ? regularized_objective
                : max_regularized_objective_raw_;
        const double elapsed =
            std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                          scale_start)
                .count();
        const double iter_rate = elapsed > 0 ? double(i + 1) / elapsed : 0.0;
        const double eta = iter_rate > 0 && nstep > i + 1
                               ? double(nstep - i - 1) / iter_rate
                               : 0.0;
        std::fprintf(
            stderr,
            "mcpd3_progress stage=dd_solve_iter scale=%ld iter=%d "
            "total_iter=%ld max_iter=%d lower_bound=%.6lf "
            "best_lower_bound=%.6lf certified_lower_bound=%.6lf "
            "best_certified_lower_bound=%.6lf regularized_objective=%.6lf "
            "best_regularized_objective=%.6lf upper_bound=%.6lf gap=%.6lf "
            "num_disagreeing=%ld disagreement_norm_sq=%.1lf "
            "step_size=%ld effective_step_size=%ld "
            "regularization_strength=%d regularization_budget=%.6lf "
            "regularization_contribution=%.6lf "
            "regularization_anchor_sink_count=%ld "
            "regularization_active_sink_count=%ld "
            "iters_since_improvement=%d solve_loop_us=%ld "
            "lagrange_update_us=%ld elapsed_sec=%.1f eta_sec=%.1f\n",
            scale_, i, total_optimization_iterations_, nstep,
            integer_to_double(lower_bound) / scale_,
            integer_to_double(best_lower_bound) / scale_,
            integer_to_double(lower_bound) / scale_,
            integer_to_double(best_lower_bound) / scale_,
            integer_to_double(regularized_objective) / scale_,
            integer_to_double(best_regularized_objective) / scale_,
            !has_current_upper_bound_
                ? std::numeric_limits<double>::infinity()
                : integer_to_double(current_upper_bound_) / scale_,
            !has_current_upper_bound_
                ? std::numeric_limits<double>::infinity()
                : integer_to_double(current_upper_bound_ - lower_bound) /
                      scale_,
            static_cast<long>(disagreeing_global_indices_.size()),
            update_stats.disagreement_norm_sq, step_size,
            update_stats.effective_step_size, regularization_strength,
            integer_to_double(last_regularization_budget_) / scale_,
            integer_to_double(last_regularization_contribution_) / scale_,
            last_regularization_anchor_sink_count_,
            last_regularization_active_sink_count_,
            i - last_improvement_iter, solve_loop_time.count(),
            lagrange_update_time.count(), elapsed, eta);
        std::fflush(stderr);
      }
      if (options_.verbose) {
        printf("iter : %6d lower_bound : %8.6lf best_lower_bound : %8.6lf upper_bound : %8.6lf gap : %8.6lf num_disagreeing : %6ld disagreement_norm_sq : %8.1lf step_size : %8ld regularization_strength : %6d regularization_budget : %8.6lf regularization_contribution : %8.6lf regularization_anchor_sink_count : %6ld regularization_active_sink_count : %6ld iters_since_improvement : %6d solve_loop_time: %8ldms lagrange_update_time: %8ldms\n",
               i, integer_to_double(lower_bound) / scale_,
               integer_to_double(best_lower_bound) / scale_,
               !has_current_upper_bound_
                   ? std::numeric_limits<double>::infinity()
                   : integer_to_double(current_upper_bound_) / scale_,
               !has_current_upper_bound_
                   ? std::numeric_limits<double>::infinity()
                   : integer_to_double(current_upper_bound_ - lower_bound) /
                         scale_,
               disagreeing_global_indices_.size(),
               update_stats.disagreement_norm_sq, update_stats.effective_step_size,
               regularization_strength,
               integer_to_double(last_regularization_budget_) / scale_,
               integer_to_double(last_regularization_contribution_) / scale_,
               last_regularization_anchor_sink_count_,
               last_regularization_active_sink_count_,
               i - last_improvement_iter, solve_loop_time.count(),
               lagrange_update_time.count());
      }

      max_lower_bound_ = std::max<double>(
          max_lower_bound_, integer_to_double(lower_bound) / scale_);
      if (!has_max_lower_bound_raw_ || lower_bound > max_lower_bound_raw_) {
        max_lower_bound_raw_ = lower_bound;
        has_max_lower_bound_raw_ = true;
      }
      if (!has_max_regularized_objective_raw_ ||
          regularized_objective > max_regularized_objective_raw_) {
        max_regularized_objective_raw_ = regularized_objective;
        has_max_regularized_objective_raw_ = true;
      }

      if (options_.iteration_callback) {
        options_.iteration_callback(DualDecompositionIterationRecord{
            total_optimization_iterations_,
            i + 1,
            scale_,
            step_size,
            update_stats.effective_step_size,
            lower_bound,
            max_lower_bound_raw_,
            regularized_objective,
            update_stats.disagreement_count,
            update_stats.disagreement_norm_sq,
            regularization_strength,
            last_regularization_budget_,
            last_regularization_contribution_,
            solve_loop_time.count(),
            lagrange_update_time.count()});
      }

      if (!has_scale_max_lower_bound || lower_bound > max_lower_bound) {
        max_lower_bound = lower_bound;
        has_scale_max_lower_bound = true;
        if (standard_early_exit_enabled &&
            !shouldSuppressEarlyScaleExit(regularization_strength) &&
            options_.legacy_patience &&
            i - last_improvement_iter >= options_.patience) {
          if (report_progress) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop reason=legacy_"
                         "patience iter=%d patience=%d\n",
                         i, options_.patience);
            std::fflush(stderr);
          }
          if (options_.verbose) {
            printf("breaking because >= %d iters since last max\n",
                   options_.patience);
          }
          opt_status = NO_FURTHER_PROGRESS;
          break;
        }
        last_improvement_iter = i;
      } else if (standard_early_exit_enabled &&
                 !shouldSuppressEarlyScaleExit(regularization_strength) &&
                 !options_.legacy_patience &&
                 i - last_improvement_iter >= options_.patience) {
        if (report_progress) {
          std::fprintf(stderr,
                       "mcpd3_progress stage=dd_solve_stop "
                       "reason=no_lower_bound_improvement iter=%d "
                       "patience=%d\n",
                       i, options_.patience);
          std::fflush(stderr);
        }
        if (options_.verbose) {
          printf("breaking because no lower-bound improvement for >= %d iters\n",
                 options_.patience);
        }
        opt_status = NO_FURTHER_PROGRESS;
        break;
      }

      if (has_best_upper_bound_ &&
          max_lower_bound >= best_upper_bound_) {
        if (report_progress) {
          if (regularization_strength == 0) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=lower_bound_closed_upper lower=%.6lf "
                         "upper=%.6lf\n",
                         integer_to_double(max_lower_bound) / scale_,
                         integer_to_double(best_upper_bound_) / scale_);
          } else {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=regularized_closed_upper lower=%.6lf "
                         "upper=%.6lf regularization_strength=%d\n",
                         integer_to_double(max_lower_bound) / scale_,
                         integer_to_double(best_upper_bound_) / scale_,
                         regularization_strength);
          }
          std::fflush(stderr);
        }
        if (options_.verbose) {
          if (regularization_strength == 0) {
            printf("breaking because lower bound closed primal upper bound: lower=%8.6lf upper=%8.6lf\n",
                   integer_to_double(max_lower_bound) / scale_,
                   integer_to_double(best_upper_bound_) / scale_);
          } else {
            printf("breaking because scaled epsilon regularization closed primal upper bound: lower=%8.6lf upper=%8.6lf regularization_strength=%d\n",
                   integer_to_double(max_lower_bound) / scale_,
                   integer_to_double(best_upper_bound_) / scale_,
                   regularization_strength);
          }
        }
        opt_status = OPTIMAL;
        break;
      }

      if (standard_early_exit_enabled && !plateau_activated_now) {
        lower_bound_group_stats.addValue(integer_to_double(lower_bound));
      }
      if (standard_early_exit_enabled && !plateau_activated_now &&
          !shouldSuppressEarlyScaleExit(regularization_strength) &&
          options_.enable_group_stopping &&
          lower_bound_group_stats.areGroupsPopulated()) {
        auto [first_group_max, second_group_max] =
            lower_bound_group_stats.getMaximums();
        if (second_group_max <= first_group_max) {
          if (report_progress) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=group_stopping first_group_max=%.6lf "
                         "second_group_max=%.6lf\n",
                         double(first_group_max) / scale_,
                         double(second_group_max) / scale_);
            std::fflush(stderr);
          }
          if (options_.verbose) {
            printf("breaking because max of this group's interval is less "
                   "than or qual to max in last last group's interval\n");
          }
          opt_status = NO_FURTHER_PROGRESS;
          break;
        }
      }

      if (disagreeing_global_indices_.size() == 0) { // optimality condition
        if (regularization_strength == 0) {
          if (report_progress) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=no_disagreement iter=%d lower=%.6lf\n",
                         i, integer_to_double(lower_bound) / scale_);
            std::fflush(stderr);
          }
          if (options_.verbose) {
            printf("breaking because of no disagreement\n");
          }
          opt_status = OPTIMAL;
        } else {
          if (report_progress) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=regularized_no_disagreement "
                         "regularization_strength=%d iter=%d lower=%.6lf\n",
                         regularization_strength, i,
                         integer_to_double(lower_bound) / scale_);
            std::fflush(stderr);
          }
          if (options_.verbose) {
            std::cout
                << "breaking because scaled epsilon regularized subproblems "
                   "agree: regularization_strength="
                << regularization_strength << " regularization_budget="
                << integer_to_string(last_regularization_budget_)
                << " regularization_budget_limit="
                << integer_to_string(regularizationBudgetLimit()) << "\n";
          }
          opt_status = OPTIMAL;
        }
        break;
      }

      //dual_cycle_list.addNode(
      //    getDualSolutionHash(disagreeing_global_indices_, lower_bound));
      //if (dual_cycle_list.getMaxCycleCount() >
      //    max_cycle_count) { // at least one set of configurations likely
      //                       // repeated more than a specified number of times
      //  printf("breaking because cycle detected\n");
      //  opt_status = NO_FURTHER_PROGRESS;
      //  break;
      //}
      if (disagreement_plateau_mode &&
          disagreement_plateau_tracker->active() &&
          !plateau_regularization_pulse_now &&
          regularization_strength != 0) {
        regularization_strength = 0;
        for (auto &solver_uptr : solvers_) {
          solver_uptr->setRegularizationStrength(regularization_strength);
        }
      }
    }
    if (options_.verbose) {
      std::cout << " === MAX === lower_bound : "
                << integer_to_string(max_lower_bound) << "\n";
    }
    if (report_progress) {
      const double elapsed =
          std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                        scale_start)
              .count();
      std::fprintf(stderr,
                   "mcpd3_progress stage=dd_solve_scale_done scale=%ld "
                   "status=%d best_lower_bound=%.6lf elapsed_sec=%.1f\n",
                   scale_, static_cast<int>(opt_status),
                   integer_to_double(max_lower_bound) / scale_, elapsed);
      std::fflush(stderr);
    }
    return opt_status;
  }

  template <int scale> void scaleProblem() {
    scaleProblem(static_cast<long>(scale));
  }

  void scaleProblem(long scale) {
    requireConstructedSolvers("scaleProblem");
    if (scale <= 0) {
      throw std::runtime_error("problem scale factor must be positive");
    }
    scale_ = checkedScaleLong(scale_, scale);
    options_.objective_scale = checkedScaleLong(options_.objective_scale, scale);
    for (auto &cap : original_arc_capacities_) {
      cap = checked_scale_capacity(cap, scale,
                                   options_.saturate_capacity_overflow);
    }
    for (auto &cap : original_terminal_capacities_) {
      cap = checked_scale_capacity(cap, scale,
                                   options_.saturate_capacity_overflow);
    }
    for (auto &solver_uptr : solvers_) {
      auto *solver = solver_uptr.get();
      const bool saturate_capacity_overflow =
          options_.saturate_capacity_overflow;
      thread_pool_.push([solver, scale, saturate_capacity_overflow] {
        solver->scaleProblem(scale, saturate_capacity_overflow);
      });
    }
    thread_pool_.wait();
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      for (auto &constraint : constraints) {
        constraint.alpha = checked_scale(
            constraint.alpha, scale, "lagrange scale promotion overflow");
        constraint.last_alpha =
            checked_scale(constraint.last_alpha, scale,
                          "lagrange scale promotion overflow");
      }
    }
    if (has_max_lower_bound_raw_) {
      max_lower_bound_raw_ = checked_scale(max_lower_bound_raw_, scale);
    }
    max_regularized_objective_raw_ = 0;
    has_max_regularized_objective_raw_ = false;
    if (has_best_upper_bound_) {
      best_upper_bound_ = checked_scale(best_upper_bound_, scale);
    }
    if (has_current_upper_bound_) {
      current_upper_bound_ = checked_scale(current_upper_bound_, scale);
    }
    warned_regularization_budget_exceeded_ = false;
  }

private:
  static long checkedScaleLong(long value, long scale) {
    if (scale <= 0) {
      throw std::runtime_error("scale factor must be positive");
    }
    if (value > 0 && value > std::numeric_limits<long>::max() / scale) {
      throw std::overflow_error("objective scale promotion overflow");
    }
    if (value < 0 && value < std::numeric_limits<long>::min() / scale) {
      throw std::overflow_error("objective scale promotion overflow");
    }
    return value * scale;
  }

  static size_t resolveThreadCount(int npartition, size_t requested) {
    size_t hardware = std::thread::hardware_concurrency();
    if (hardware == 0) {
      hardware = 1;
    }
    size_t limit = requested == 0 ? hardware : requested;
    return std::max<size_t>(1, std::min<size_t>(npartition, limit));
  }

  long getDualSolutionHash(const std::list<int> &disagreeing_global_indices,
                           const Objective &lower_bound) const {
    std::hash<std::string> hasher{};
    long hash = static_cast<long>(hasher(integer_to_string(lower_bound)));
    for (const auto &global_index : disagreeing_global_indices) {
      hash ^= static_cast<long>(hasher(std::to_string(global_index)));
    }
    return hash;
  }

  LagrangeUpdateStats runLagrangeMultipliersUpdateStep(long step_size,
                                                       bool use_momentum,
                                                       const Objective &lower_bound) {
    LagrangeUpdateStats stats;
    for (auto &[global_index, constraints] : constraint_arc_map_) {
      bool disagreement_exists = false;
      for (auto &constraint : constraints) {
        int diff =
            solvers_[constraint.partition_index_target]->getMinCutSolution(
                constraint.local_index_target) -
            solvers_[constraint.partition_index_source]->getMinCutSolution(
                constraint.local_index_source);
        if (diff != 0) {
          disagreement_exists = true;
          stats.disagreement_count += std::abs(diff);
          stats.disagreement_norm_sq += static_cast<double>(diff * diff);
        }
      }
      if (disagreement_exists) {
        stats.disagreeing_global_indices.emplace_back(global_index);
      }
    }
    (void)lower_bound;
    stats.effective_step_size = std::max(step_size, options_.min_step_size);

    for (auto &[global_index, constraints] : constraint_arc_map_) {
      for (auto &constraint : constraints) {
        constraint.last_alpha = constraint.alpha; // record alpha before update
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
            const long alpha_update =
                stats.effective_step_size *
                static_cast<int>(momentum_scale * constraint.alpha_momentum);
            constraint.alpha = checked_add(
                constraint.alpha, lagrange_from_integer(alpha_update),
                "lagrange multiplier overflow");
          } else {
            constraint.alpha = checked_add(
                constraint.alpha,
                lagrange_from_integer(stats.effective_step_size * diff),
                "lagrange multiplier overflow");
          }
        }
      }
    }
    return stats;
  }

  void validateOptions() const {
    if (options_.objective_scale <= 0) {
      throw std::runtime_error("objective scale must be positive");
    }
    if (options_.halo_depth != kInfiniteHaloDepth &&
        options_.halo_depth < 1) {
      throw std::runtime_error(
          "halo depth must be positive or kInfiniteHaloDepth");
    }
    if (options_.initial_alpha_random_radius < 0) {
      throw std::runtime_error(
          "initial alpha random radius must be non-negative");
    }
    if (options_.regularization_budget_limit < 0) {
      throw std::runtime_error(
          "regularization budget limit must be non-negative");
    }
    if (options_.scaled_epsilon_max_step_size <= 0) {
      throw std::runtime_error(
          "scaled epsilon maximum step size must be positive");
    }
    if (options_.scaled_epsilon_strength_cap < 0) {
      throw std::runtime_error(
          "scaled epsilon strength cap must be non-negative");
    }
    if (options_.regularization_scheme ==
            DualDecompositionRegularizationScheme::
                DISAGREEMENT_PLATEAU_EPSILON &&
        options_.disagreement_patience <= 0) {
      throw std::runtime_error("disagreement patience must be positive");
    }
    if (options_.max_objective_scale_promotions < 0) {
      throw std::runtime_error(
          "max objective scale promotions must be non-negative");
    }
    if (options_.max_total_iteration_count < 0) {
      throw std::runtime_error(
          "maximum total iteration count must be non-negative");
    }
    if (!options_.partition_labels.empty()) {
      if (!options_.partition_edge_weights.empty()) {
        throw std::runtime_error(
            "partition labels and partition edge weights are mutually "
            "exclusive");
      }
      if (options_.partition_labels.size() !=
          static_cast<size_t>(nnode_)) {
        throw std::runtime_error(
            "partition label count must match the global node count");
      }
      std::vector<bool> populated(static_cast<size_t>(npartition_), false);
      for (const int label : options_.partition_labels) {
        if (label < 0 || label >= npartition_) {
          throw std::runtime_error(
              "partition label must be within the partition range");
        }
        populated[static_cast<size_t>(label)] = true;
      }
      if (std::find(populated.begin(), populated.end(), false) !=
          populated.end()) {
        throw std::runtime_error(
            "partition labels must populate every partition");
      }
    }
    if (!options_.reference_cut_labels.empty()) {
      if (options_.reference_cut_check_interval <= 0) {
        throw std::runtime_error(
            "reference cut check interval must be positive");
      }
      if (options_.reference_cut_labels.size() !=
          static_cast<size_t>(nnode_)) {
        throw std::runtime_error(
            "reference cut label count must match the global node count");
      }
      if (options_.canonical_cut_selection !=
          CanonicalCutSelection::SOLVER_DEFAULT) {
        throw std::runtime_error(
            "reference and canonical cut selection are mutually exclusive");
      }
      for (const int label : options_.reference_cut_labels) {
        if (label != 0 && label != 1) {
          throw std::runtime_error("reference cut labels must be binary");
        }
      }
    }
    if (!options_.construct_solvers) {
      if (!options_.emit_partition_packages) {
        throw std::runtime_error(
            "solver construction can only be disabled when partition package "
            "export is enabled");
      }
      if (options_.track_primal_upper_bound) {
        throw std::runtime_error(
            "solver construction can only be disabled when primal upper bound "
            "tracking is disabled");
      }
    }
  }

  void requireConstructedSolvers(const char *operation) const {
    if (!options_.construct_solvers) {
      throw std::runtime_error(std::string(operation) +
                               " requires constructed solvers");
    }
  }

  Objective regularizationBudgetLimit() const {
    return options_.regularization_budget_limit > 0
               ? options_.regularization_budget_limit
               : Objective(options_.objective_scale);
  }

  bool totalIterationBudgetExhausted() const {
    return options_.max_total_iteration_count > 0 &&
           total_optimization_iterations_ >=
               options_.max_total_iteration_count;
  }

  void warnIfRegularizationBudgetExceeded(const Objective &budget,
                                          const Capacity &regularization_strength) {
    if (!isRegularizationBudgetExceeded(budget, regularization_strength) ||
        warned_regularization_budget_exceeded_) {
      return;
    }
    std::fprintf(
        stderr,
        "warning: regularization budget %s is not below limit %s; "
        "a regularized agreement may not certify optimality\n",
        integer_to_string(budget).c_str(),
        integer_to_string(regularizationBudgetLimit()).c_str());
    std::fflush(stderr);
    warned_regularization_budget_exceeded_ = true;
  }

  bool isRegularizationBudgetExceeded(
      const Objective &budget, const Capacity &regularization_strength) const {
    (void)regularization_strength;
    return budget >= regularizationBudgetLimit();
  }

  bool shouldSuppressEarlyScaleExit(
      const Capacity &regularization_strength) const {
    return options_.exhaust_scale_iterations ||
           (options_.exhaust_regularized_scale_iterations &&
            regularization_strength > 0);
  }

  bool tryPromoteObjectiveScale(long factor, long *step_size) {
    if (!options_.promote_objective_scale_on_overbudget ||
        options_.regularization_budget_limit > 0 ||
        objective_scale_promotion_count_ >=
            options_.max_objective_scale_promotions) {
      return false;
    }
    const long old_scale = scale_;
    scaleProblem(factor);
    ++objective_scale_promotion_count_;
    *step_size = scale_;
    options_.initial_step_size = *step_size;
    if (dualdecomp_progress_enabled()) {
      std::fprintf(stderr,
                   "mcpd3_progress stage=dd_objective_scale_promote "
                   "old_scale=%ld new_scale=%ld factor=%ld "
                   "restart_step_size=%ld promotion_count=%ld\n",
                   old_scale, scale_, factor, *step_size,
                   objective_scale_promotion_count_);
      std::fflush(stderr);
    }
    if (options_.verbose) {
      printf("promoting objective scale from %ld to %ld and restarting at step %ld\n",
             old_scale, scale_, *step_size);
    }
    return true;
  }

  Objective computePrimalCutValue(const std::vector<bool> &labels) const {
    Objective cut_value = 0;
    for (int i = 0; i < narc_; ++i) {
      const int s = original_arcs_[2 * i + 0];
      const int t = original_arcs_[2 * i + 1];
      const Capacity forward_capacity = original_arc_capacities_[2 * i + 0];
      const Capacity backward_capacity = original_arc_capacities_[2 * i + 1];
      if (!labels[s] && labels[t]) {
        cut_value = checked_add(
            cut_value, widen_capacity(forward_capacity),
            "primal cut objective overflow");
      } else if (labels[s] && !labels[t]) {
        cut_value = checked_add(
            cut_value, widen_capacity(backward_capacity),
            "primal cut objective overflow");
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      const Capacity terminal_capacity = original_terminal_capacities_[i];
      if (!labels[i] && terminal_capacity < 0) {
        cut_value = checked_add(
            cut_value, absolute_capacity(terminal_capacity),
            "primal terminal objective overflow");
      } else if (labels[i] && terminal_capacity > 0) {
        cut_value = checked_add(
            cut_value, widen_capacity(terminal_capacity),
            "primal terminal objective overflow");
      }
    }
    return cut_value;
  }

  Objective updatePrimalUpperBound() {
    std::vector<int> vote_count(nnode_, 0);
    std::vector<int> sink_vote_count(nnode_, 0);
    for (int i = 0; i < npartition_; ++i) {
      const auto &min_cut_sub_graph = min_cut_sub_graphs_[i];
      const auto &solver = solvers_[i];
      for (int local_index = 0;
           local_index < static_cast<int>(min_cut_sub_graph.local_to_global.size());
           ++local_index) {
        const int global_index = min_cut_sub_graph.local_to_global[local_index];
        ++vote_count[global_index];
        sink_vote_count[global_index] +=
            solver->getMinCutSolution(local_index) ? 1 : 0;
      }
    }
    std::vector<bool> decoded(nnode_, false);
    for (int i = 0; i < nnode_; ++i) {
      // Deterministic tie-break: source side, label 0.
      decoded[i] = sink_vote_count[i] * 2 > vote_count[i];
    }
    const Objective upper_bound = computePrimalCutValue(decoded);
    if (!has_best_upper_bound_ || upper_bound < best_upper_bound_) {
      best_upper_bound_ = upper_bound;
      has_best_upper_bound_ = true;
      best_primal_solution_ = decoded;
    }
    primal_solution_ = decoded;
    return upper_bound;
  }

  void validateAndReportPartition(const std::vector<int> &partitions) const {
    if (static_cast<int>(partitions.size()) != nnode_) {
      throw std::runtime_error("partition vector size does not match node count");
    }

    std::vector<long> part_node_counts(npartition_, 0);
    std::vector<long> part_arc_counts(npartition_, 0);
    std::vector<unsigned char> boundary_nodes(nnode_, 0);
    std::vector<unsigned char> constrained_nodes(nnode_, 0);

    for (int node = 0; node < nnode_; ++node) {
      const int part = partitions[node];
      if (part < 0 || part >= npartition_) {
        throw std::runtime_error("partition label out of range");
      }
      ++part_node_counts[part];
    }

    long crossing_edges = 0;
    long boundary_node_count = 0;
    long constrained_node_count = 0;
    for (int aid = 0; aid < narc_; ++aid) {
      int s = arcs_[2 * aid + 0];
      int t = arcs_[2 * aid + 1];
      if (s < 0 || s >= nnode_ || t < 0 || t >= nnode_) {
        throw std::runtime_error("arc endpoint out of range");
      }
      if (s > t) {
        std::swap(s, t);
      }
      const int s_part = partitions[s];
      const int t_part = partitions[t];
      ++part_arc_counts[s_part];
      if (s_part != t_part) {
        ++crossing_edges;
        if (!boundary_nodes[s]) {
          boundary_nodes[s] = 1;
          ++boundary_node_count;
        }
        if (!boundary_nodes[t]) {
          boundary_nodes[t] = 1;
          ++boundary_node_count;
        }
        // This mirrors initializeDecomposition(): after endpoint ordering, the
        // arc belongs to s_part and t becomes the constrained clone node.
        if (!constrained_nodes[t]) {
          constrained_nodes[t] = 1;
          ++constrained_node_count;
        }
      }
    }

    const auto [min_nodes_it, max_nodes_it] =
        std::minmax_element(part_node_counts.begin(), part_node_counts.end());
    const auto [min_arcs_it, max_arcs_it] =
        std::minmax_element(part_arc_counts.begin(), part_arc_counts.end());
    const double mean_nodes =
        std::accumulate(part_node_counts.begin(), part_node_counts.end(), 0.0) /
        static_cast<double>(std::max(1, npartition_));
    const double mean_arcs =
        std::accumulate(part_arc_counts.begin(), part_arc_counts.end(), 0.0) /
        static_cast<double>(std::max(1, npartition_));
    const double node_imbalance =
        mean_nodes > 0.0 ? static_cast<double>(*max_nodes_it) / mean_nodes
                         : 0.0;
    const double arc_imbalance =
        mean_arcs > 0.0 ? static_cast<double>(*max_arcs_it) / mean_arcs : 0.0;

    if (dualdecomp_progress_enabled()) {
      std::fprintf(
          stderr,
          "mcpd3_progress stage=partition_validation nnode=%d narc=%d "
          "partitions=%d crossing_edges=%ld boundary_nodes=%ld "
          "constrained_nodes=%ld min_part_nodes=%ld max_part_nodes=%ld "
          "mean_part_nodes=%.1f node_imbalance=%.3f min_part_arcs=%ld "
          "max_part_arcs=%ld mean_part_arcs=%.1f arc_imbalance=%.3f\n",
          nnode_, narc_, npartition_, crossing_edges, boundary_node_count,
          constrained_node_count, *min_nodes_it, *max_nodes_it, mean_nodes,
          node_imbalance, *min_arcs_it, *max_arcs_it, mean_arcs,
          arc_imbalance);
      std::fflush(stderr);
    }
  }

  void initializeDecomposition() {
    const bool report_progress = dualdecomp_progress_enabled();
    const long progress_interval = 10000000;
    const auto init_start = std::chrono::steady_clock::now();
    std::unordered_map</*global_index=*/int,
                       /*exists_in_partitions=*/std::set<int>>
        constrained_nodes;
    /**
     * step 0: parition graph into npartition_ partitions
     */
    const auto *partition_edge_weights =
        options_.partition_edge_weights.empty()
            ? nullptr
            : &options_.partition_edge_weights;
    partition_labels_ = options_.partition_labels.empty()
                            ? configured_graph_partition(
                                  npartition_, narc_, nnode_, arcs_,
                                  &arc_capacities_, partition_edge_weights)
                            : options_.partition_labels;
    const auto &partitions_ = partition_labels_;
    validateAndReportPartition(partitions_);
    HaloPartitionLayout halo_layout;
    if (options_.halo_depth != 1) {
      halo_layout = buildHaloPartitionLayout(
          npartition_, nnode_, arcs_, partition_labels_, options_.halo_depth);
    }
    halo_objective_multiplier_ = halo_layout.objective_multiplier;
    if (halo_objective_multiplier_ > 1) {
      options_.objective_scale = checkedScaleLong(
          options_.objective_scale, halo_objective_multiplier_);
      if (options_.regularization_budget_limit > 0) {
        options_.regularization_budget_limit = checked_scale(
            options_.regularization_budget_limit, halo_objective_multiplier_,
            "halo regularization budget overflow");
      }
      for (Capacity &capacity : original_arc_capacities_) {
        capacity = checked_scale_capacity(
            capacity, halo_objective_multiplier_,
            options_.saturate_capacity_overflow);
      }
      for (Capacity &capacity : original_terminal_capacities_) {
        capacity = checked_scale_capacity(
            capacity, halo_objective_multiplier_,
            options_.saturate_capacity_overflow);
      }
    }
    dualdecomp_progress_report("dd_partition_done", 1, 1, init_start);
    auto mapping_start = std::chrono::steady_clock::now();
    int mapping_done = 0;
    for (auto &min_cut_sub_graph : min_cut_sub_graphs_) {
      min_cut_sub_graph.initializeMapping(nnode_);
      ++mapping_done;
      dualdecomp_progress_report("dd_initialize_mapping", mapping_done,
                                 npartition_, mapping_start);
    }
    /** step 1: distribute arcs into their local halo subgraphs. */
    auto arc_start = std::chrono::steady_clock::now();
    if (options_.halo_depth == 1) {
      arc_locations_.resize(static_cast<size_t>(narc_));
      for (int i = 0; i < narc_; ++i) {
        int s = arcs_[2 * i + 0];
        int t = arcs_[2 * i + 1];
        Capacity forward_capacity = arc_capacities_[2 * i + 0];
        Capacity backward_capacity = arc_capacities_[2 * i + 1];
        bool swapped = false;
        if (s > t) {
          std::swap(s, t);
          std::swap(forward_capacity, backward_capacity);
          swapped = true;
        }
        const int arc_partition = partitions_[s];
        auto &min_cut_sub_graph = min_cut_sub_graphs_[arc_partition];
        arc_locations_[static_cast<size_t>(i)] = ArcLocation{
            arc_partition, min_cut_sub_graph.graph.narc, swapped};
        min_cut_sub_graph.insertArc(s, t, forward_capacity,
                                    backward_capacity);
        if (partitions_[t] != arc_partition) {
          constrained_nodes[t].insert(arc_partition);
        }
        if (report_progress && (i + 1) % progress_interval == 0) {
          dualdecomp_progress_report("dd_distribute_arcs", i + 1, narc_,
                                     arc_start);
        }
      }
    } else {
      halo_arc_locations_.resize(static_cast<size_t>(narc_));
      halo_arc_objective_factors_.resize(static_cast<size_t>(narc_));
      for (int i = 0; i < narc_; ++i) {
        int s = arcs_[2 * i + 0];
        int t = arcs_[2 * i + 1];
        Capacity forward_capacity = arc_capacities_[2 * i + 0];
        Capacity backward_capacity = arc_capacities_[2 * i + 1];
        bool swapped = false;
        if (s > t) {
          std::swap(s, t);
          std::swap(forward_capacity, backward_capacity);
          swapped = true;
        }
        const auto &arc_partitions =
            halo_layout.arc_partitions[static_cast<size_t>(i)];
        const long objective_factor =
            halo_objective_multiplier_ /
            static_cast<long>(arc_partitions.size());
        halo_arc_objective_factors_[static_cast<size_t>(i)] = objective_factor;
        for (const int arc_partition : arc_partitions) {
          auto &min_cut_sub_graph = min_cut_sub_graphs_[arc_partition];
          halo_arc_locations_[static_cast<size_t>(i)].push_back(ArcLocation{
              arc_partition, min_cut_sub_graph.graph.narc, swapped});
          min_cut_sub_graph.insertArc(
              s, t,
              checked_scale_capacity(forward_capacity, objective_factor,
                                     options_.saturate_capacity_overflow),
              checked_scale_capacity(backward_capacity, objective_factor,
                                     options_.saturate_capacity_overflow));
        }
        if (report_progress && (i + 1) % progress_interval == 0) {
          dualdecomp_progress_report("dd_distribute_arcs", i + 1, narc_,
                                     arc_start);
        }
      }
    }
    dualdecomp_progress_report("dd_distribute_arcs", narc_, narc_, arc_start);
    arcs_.clear();
    arcs_.shrink_to_fit();
    arc_capacities_.clear();
    arc_capacities_.shrink_to_fit();
    /** step 2: distribute node unaries and materialize halo copies. */
    auto terminal_start = std::chrono::steady_clock::now();
    if (options_.halo_depth == 1) {
      terminal_locations_.resize(static_cast<size_t>(nnode_));
      if (options_.materialize_all_partition_nodes) {
        for (int node = 0; node < nnode_; ++node) {
          min_cut_sub_graphs_[partitions_[node]].getOrInsertNode(node);
        }
      }
      for (int i = 0; i < nnode_; ++i) {
        if (terminal_capacities_[i] != 0) {
          min_cut_sub_graphs_[partitions_[i]].insertTerminal(
              i, terminal_capacities_[i]);
        }
        if (report_progress && (i + 1) % progress_interval == 0) {
          dualdecomp_progress_report("dd_distribute_terminals", i + 1,
                                     nnode_, terminal_start);
        }
      }
      for (auto &[global_index, memberships] : constrained_nodes) {
        memberships.insert(partitions_[global_index]);
        for (const int partition : memberships) {
          min_cut_sub_graphs_[partition].getOrInsertNode(global_index);
        }
      }
      for (int node = 0; node < nnode_; ++node) {
        const int partition = partitions_[node];
        const auto &global_to_local =
            min_cut_sub_graphs_[partition].global_to_local_map;
        if (node < static_cast<int>(global_to_local.size()) &&
            global_to_local[node] >= 0) {
          terminal_locations_[static_cast<size_t>(node)] =
              TerminalLocation{partition, global_to_local[node]};
        }
      }
    } else {
      halo_terminal_locations_.resize(static_cast<size_t>(nnode_));
      halo_terminal_objective_factors_.resize(static_cast<size_t>(nnode_));
      for (int node = 0; node < nnode_; ++node) {
        const auto &memberships =
            halo_layout.node_partitions[static_cast<size_t>(node)];
        const long objective_factor =
            halo_objective_multiplier_ /
            static_cast<long>(memberships.size());
        halo_terminal_objective_factors_[static_cast<size_t>(node)] =
            objective_factor;
        if (memberships.size() > 1) {
          constrained_nodes[node].insert(memberships.begin(),
                                         memberships.end());
        }
        for (const int partition : memberships) {
          const int local_node =
              min_cut_sub_graphs_[partition].getOrInsertNode(node);
          halo_terminal_locations_[static_cast<size_t>(node)].push_back(
              TerminalLocation{partition, local_node});
          if (terminal_capacities_[node] != 0) {
            min_cut_sub_graphs_[partition].insertTerminal(
                node, checked_scale_capacity(
                          terminal_capacities_[node], objective_factor,
                          options_.saturate_capacity_overflow));
          }
        }
        if (report_progress && (node + 1) % progress_interval == 0) {
          dualdecomp_progress_report("dd_distribute_terminals", node + 1,
                                     nnode_, terminal_start);
        }
      }
    }
    dualdecomp_progress_report("dd_distribute_terminals", nnode_, nnode_,
                               terminal_start);
    auto finalize_start = std::chrono::steady_clock::now();
    int finalize_done = 0;
    for (auto &min_cut_sub_graph : min_cut_sub_graphs_) {
      min_cut_sub_graph.finalizeTerminals();
      ++finalize_done;
      dualdecomp_progress_report("dd_finalize_terminals", finalize_done,
                                 npartition_, finalize_start);
    }
    local_arc_counts_.reserve(min_cut_sub_graphs_.size());
    local_node_counts_.reserve(min_cut_sub_graphs_.size());
    for (const auto &min_cut_sub_graph : min_cut_sub_graphs_) {
      local_arc_counts_.push_back(min_cut_sub_graph.graph.narc);
      local_node_counts_.push_back(min_cut_sub_graph.graph.nnode);
    }
    terminal_capacities_.clear();
    terminal_capacities_.shrink_to_fit();
    /**
     * step 3: create solvers
     */
    auto solver_start = std::chrono::steady_clock::now();
    int solver_done = 0;
    for (int partition = 0; partition < npartition_; ++partition) {
      auto &min_cut_sub_graph = min_cut_sub_graphs_[partition];
      if (options_.emit_partition_packages) {
        auto &package = partition_packages_[partition];
        package.partition_id = partition;
        package.local_node_count = min_cut_sub_graph.graph.nnode;
        package.objective_multiplier = halo_objective_multiplier_;
        if (options_.construct_solvers) {
          package.arcs = min_cut_sub_graph.graph.arcs;
          package.arc_capacities = min_cut_sub_graph.graph.arc_capacities;
          package.terminal_capacities =
              min_cut_sub_graph.graph.terminal_capacities;
          package.local_to_global = min_cut_sub_graph.local_to_global;
        } else {
          package.arcs = std::move(min_cut_sub_graph.graph.arcs);
          package.arc_capacities =
              std::move(min_cut_sub_graph.graph.arc_capacities);
          package.terminal_capacities =
              std::move(min_cut_sub_graph.graph.terminal_capacities);
          package.local_to_global =
              std::move(min_cut_sub_graph.local_to_global);
        }
        package.constraint_endpoints.clear();
      }

      if (options_.construct_solvers) {
        auto solver = std::make_unique<PrimalDualMinCutSolver>(
            std::move(min_cut_sub_graph.graph));
        solver->setTrackArcFlowUpdates(options_.track_arc_flow_updates);
        solver->setCanonicalCutSelection(options_.canonical_cut_selection);
        if (!options_.reference_cut_labels.empty()) {
          std::vector<int> local_reference;
          local_reference.reserve(min_cut_sub_graph.local_to_global.size());
          for (const int global_node : min_cut_sub_graph.local_to_global) {
            local_reference.push_back(options_.reference_cut_labels[
                static_cast<size_t>(global_node)]);
          }
          solver->setReferenceCutLabels(std::move(local_reference));
          solver->setReferenceCutSelection(options_.reference_cut_selection);
          solver->setReferenceCutCheckInterval(
              options_.reference_cut_check_interval);
        }
        solver->setForceFullMinCutRecompute(
            options_.force_full_mincut_recompute);
        solvers_.emplace_back(std::move(solver));
      }
      ++solver_done;
      dualdecomp_progress_report("dd_create_solvers", solver_done, npartition_,
                                 solver_start);
    }
    /**
     * step 4: create a DualDecompositionConstraintArc for each constraint
     * induced on each constrained node
     */
    std::set<int> constrained_nodes_partition_counts;
    std::vector<int> constrained_nodes_count_in_each_partition(npartition_,0);
    int next_constraint_id = 0;
    auto constraint_start = std::chrono::steady_clock::now();
    long constrained_done = 0;
    const long constrained_total = static_cast<long>(constrained_nodes.size());
    std::mt19937 initial_alpha_generator(options_.initial_alpha_random_seed);
    std::uniform_int_distribution<long> initial_alpha_distribution(
        -options_.initial_alpha_random_radius,
        options_.initial_alpha_random_radius);
    auto initial_alpha = [&]() -> Lagrange {
      if (!options_.randomize_initial_alphas ||
          options_.initial_alpha_random_radius == 0) {
        return 0;
      }
      return lagrange_from_integer(
          initial_alpha_distribution(initial_alpha_generator));
    };
    for (auto &[global_index, partitions] : constrained_nodes) {
      constraint_arc_map_.push_back({global_index, {}});
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
          const int constraint_id = next_constraint_id++;
          const Lagrange alpha = initial_alpha();
          constraint_arcs.emplace_back(
              /*alpha=*/alpha,
              /*last_alpha=*/alpha,
              /*alpha_momentum=*/0,
              /*partition_index_source=*/partition_source,
              /*partition_index_target=*/partition_target,
              /*local_index_source=*/local_index_source,
              /*local_index_target=*/local_index_target);
          auto arc_reference = --constraint_arcs.end();
          if (options_.construct_solvers) {
            solvers_[partition_source]->addSourceDualDecompositionConstraint(
                arc_reference);
            solvers_[partition_target]->addTargetDualDecompositionConstraint(
                arc_reference);
          }
          if (options_.emit_partition_packages) {
            partition_packages_[partition_source]
                .constraint_endpoints.push_back(
                    ConstraintEndpointBinding{/*constraint_id=*/constraint_id,
                                              /*global_node_id=*/global_index,
                                              /*local_index=*/local_index_source,
                                              /*is_source=*/true,
                                              /*alpha=*/alpha,
                                              /*last_alpha=*/alpha,
                                              /*alpha_momentum=*/0});
            partition_packages_[partition_target]
                .constraint_endpoints.push_back(
                    ConstraintEndpointBinding{/*constraint_id=*/constraint_id,
                                              /*global_node_id=*/global_index,
                                              /*local_index=*/local_index_target,
                                              /*is_source=*/false,
                                              /*alpha=*/alpha,
                                              /*last_alpha=*/alpha,
                                              /*alpha_momentum=*/0});
          }
          constrained_nodes_count_in_each_partition[partition_source]++;
          constrained_nodes_count_in_each_partition[partition_target]++;
        }
      }
      constrained_nodes_partition_counts.insert(partitions.size());
      ++constrained_done;
      if (report_progress && constrained_done % 1000000 == 0) {
        dualdecomp_progress_report("dd_create_constraints", constrained_done,
                                   constrained_total, constraint_start);
      }
    }
    dualdecomp_progress_report("dd_create_constraints", constrained_done,
                               constrained_total, constraint_start);
    constrained_nodes.clear();
    constraint_arc_map_.shrink_to_fit();
    if (!options_.track_primal_upper_bound) {
      auto release_start = std::chrono::steady_clock::now();
      int release_done = 0;
      for (auto &min_cut_sub_graph : min_cut_sub_graphs_) {
        min_cut_sub_graph.releaseConstructionMaps();
        ++release_done;
        dualdecomp_progress_report("dd_release_construction_maps",
                                   release_done, npartition_, release_start);
      }
    }
    if (!options_.track_arc_flow_updates &&
        options_.partition_edge_weights.empty()) {
      partition_labels_.clear();
      partition_labels_.shrink_to_fit();
    }
    if (!options_.construct_solvers) {
      constraint_arc_map_.clear();
      constraint_arc_map_.shrink_to_fit();
      min_cut_sub_graphs_.clear();
      min_cut_sub_graphs_.shrink_to_fit();
      primal_solution_.clear();
      primal_solution_.shrink_to_fit();
      best_primal_solution_.clear();
      best_primal_solution_.shrink_to_fit();
    }
    dualdecomp_progress_report("dd_initialize_decomposition_total", 1, 1,
                               init_start);
    if (options_.verbose) {
      printf("partition counts: ");
      for (const auto &count : constrained_nodes_partition_counts) {
        printf("%d,", count);
      }
      printf("\n");
      printf("max count of constrainted nodes in any one partition: %d\n",
          *std::max_element(constrained_nodes_count_in_each_partition.begin(),
            constrained_nodes_count_in_each_partition.end()));
      printf("mean count of constrainted nodes in any one partition: %lf\n",
             std::accumulate(constrained_nodes_count_in_each_partition.begin(),
                             constrained_nodes_count_in_each_partition.end(), 0) /
                 (static_cast<double>(
                     constrained_nodes_count_in_each_partition.size())));
    }
    /**
     * step 5: create a min cut problem from the original problem to evaluate
     * the primal objective value
     */
    // primal_solver_ = std::make_unique<PrimalDualMinCutSolver>(nnode_,narc_,
    //   std::move(arcs_),
    //   std::move(arc_capacities_),
    //   std::move(terminal_capacities_));
  }

  /**
   * data passed into decomposition
   */
  int npartition_;
  int nnode_;
  int narc_;
  std::vector<int> arcs_;
  std::vector<Capacity> arc_capacities_;
  std::vector<Capacity> terminal_capacities_;
  std::vector<int> original_arcs_;
  std::vector<Capacity> original_arc_capacities_;
  std::vector<Capacity> original_terminal_capacities_;
  std::vector<int> partition_labels_;

  struct ArcLocation {
    int partition = -1;
    int local_arc = -1;
    bool swapped = false;
  };
  struct TerminalLocation {
    int partition = -1;
    int local_node = -1;
  };
  std::vector<ArcLocation> arc_locations_;
  std::vector<TerminalLocation> terminal_locations_;
  std::vector<std::vector<ArcLocation>> halo_arc_locations_;
  std::vector<std::vector<TerminalLocation>> halo_terminal_locations_;
  std::vector<long> halo_arc_objective_factors_;
  std::vector<long> halo_terminal_objective_factors_;
  std::vector<int> local_arc_counts_;
  std::vector<int> local_node_counts_;

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
    std::vector</*local_index -> global_index*/ int> local_to_global;
    std::vector</*global_index -> local_index*/ int> global_to_local_map;

    MinCutSubGraph() {
      graph.nnode = 0;
      graph.narc = 0;
    }

    void initializeMapping(int global_node_count) {
      global_to_local_map.assign(global_node_count, -1);
      local_to_global.clear();
    }

    int getOrInsertNode(int global_index) {
      if (global_index < 0 ||
          global_index >= static_cast<int>(global_to_local_map.size())) {
        throw std::runtime_error("global node index out of range");
      }
      int &local_index = global_to_local_map[global_index];
      if (local_index < 0) {
        local_index = graph.nnode++;
        local_to_global.push_back(global_index);
      }
      return local_index;
    }

    int getNode(int global_index) const {
      if (global_index < 0 ||
          global_index >= static_cast<int>(global_to_local_map.size()) ||
          global_to_local_map[global_index] < 0) {
        throw std::runtime_error("Node not found");
      }
      return global_to_local_map[global_index];
    }

    void insertArc(int global_source_index, int global_target_index,
                   const Capacity &forward_capacity,
                   const Capacity &backward_capacity) {
      int s = getOrInsertNode(global_source_index);
      int t = getOrInsertNode(global_target_index);
      graph.arc_capacities.push_back(forward_capacity);
      graph.arc_capacities.push_back(backward_capacity);
      graph.arcs.push_back(s);
      graph.arcs.push_back(t);
      graph.narc++;
    }

    void insertTerminal(int global_index, const Capacity &terminal_capacity) {
      int s = getOrInsertNode(global_index);
      if (static_cast<int>(graph.terminal_capacities.size()) < graph.nnode) {
        graph.terminal_capacities.resize(graph.nnode, 0);
      }
      graph.terminal_capacities[s] = terminal_capacity;
    }

    void finalizeTerminals() {
      graph.terminal_capacities.resize(graph.nnode, 0);
    }

    void releaseConstructionMaps() {
      local_to_global.clear();
      local_to_global.shrink_to_fit();
      global_to_local_map.clear();
      global_to_local_map.shrink_to_fit();
    }
  };

  std::vector<MinCutSubGraph> min_cut_sub_graphs_;
  std::vector<PartitionPackage> partition_packages_;
  std::vector<bool> primal_solution_;
  std::vector<bool> best_primal_solution_;
  long scale_;
  DualDecompositionOptions options_;
  ThreadPool<void> thread_pool_;
  long solve_loop_time_;
  long lagrange_update_time_;
  double max_lower_bound_;
  Objective max_lower_bound_raw_;
  Objective max_regularized_objective_raw_;
  Objective best_upper_bound_;
  Objective current_upper_bound_;
  bool has_max_lower_bound_raw_;
  bool has_max_regularized_objective_raw_;
  bool has_best_upper_bound_;
  bool has_current_upper_bound_;
  Objective last_original_objective_raw_;
  Objective last_certified_lower_bound_raw_;
  Objective last_regularized_objective_raw_;
  long last_disagreement_count_;
  double last_disagreement_norm_sq_;
  Objective last_regularization_budget_;
  Objective last_regularization_contribution_;
  long last_regularization_anchor_sink_count_;
  long last_regularization_active_sink_count_;
  long total_optimization_iterations_;
  long unit_step_no_momentum_retry_count_ = 0;
  long objective_scale_promotion_count_;
  long halo_objective_multiplier_;
  long disagreement_plateau_activation_count_ = 0;
  bool warned_regularization_budget_exceeded_;
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
