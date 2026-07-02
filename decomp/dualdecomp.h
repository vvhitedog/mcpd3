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
#include <cstdio>
#include <cstdlib>
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
#include <decomp/lower_bound_certificate.h>
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
  NONE
};

struct DualDecompositionOptions {
  int num_optimization_scales = 5;
  int max_iteration_count = 10000;
  int max_cycle_count = 2;
  long initial_step_size = 10000;
  int patience = 10;
  bool legacy_patience = false;
  bool use_momentum = true;
  bool enable_group_stopping = true;
  bool track_primal_upper_bound = true;
  bool emit_partition_packages = true;
  bool saturate_capacity_overflow = false;
  bool verbose = true;
  long min_step_size = 1;
  long max_step_size = 10000;
  long objective_scale = 1;
  size_t thread_count = 0;
  DualDecompositionRegularizationScheme regularization_scheme =
      DualDecompositionRegularizationScheme::SCALED_EPSILON;
  long regularization_budget_limit = 0;
  bool promote_objective_scale_on_overbudget = true;
  int max_objective_scale_promotions = 4;
  bool randomize_initial_alphas = false;
  long initial_alpha_random_radius = 0;
  unsigned int initial_alpha_random_seed = 0;
};

class DualDecomposition {
public:
  DualDecomposition(int npartition, int nnode, int narc, std::vector<int> arcs,
                    std::vector<int> arc_capacities,
                    std::vector<int> terminal_capacities,
                    DualDecompositionOptions options = {})
      : npartition_(npartition), nnode_(nnode), narc_(narc),
        arcs_(std::move(arcs)), arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)),
        original_arcs_(options.track_primal_upper_bound ? arcs_
                                                        : std::vector<int>()),
        original_arc_capacities_(options.track_primal_upper_bound
                                     ? arc_capacities_
                                     : std::vector<int>()),
        original_terminal_capacities_(options.track_primal_upper_bound
                                          ? terminal_capacities_
                                          : std::vector<int>()),
        min_cut_sub_graphs_(npartition_),
        partition_packages_(options.emit_partition_packages ? npartition_ : 0),
        primal_solution_(nnode_), scale_(1),
        options_(options),
        thread_pool_(
            resolveThreadCount(npartition_, options_.thread_count)),
        solve_loop_time_(0),
        max_lower_bound_(std::numeric_limits<double>::lowest()),
        max_lower_bound_raw_(std::numeric_limits<long>::min()),
        max_regularized_objective_raw_(std::numeric_limits<long>::min()),
        best_upper_bound_(std::numeric_limits<long>::max()),
        current_upper_bound_(std::numeric_limits<long>::max()),
        last_disagreement_count_(0),
        last_disagreement_norm_sq_(0),
        last_regularization_budget_(0),
        last_regularization_contribution_(0),
        last_regularization_anchor_sink_count_(0),
        last_regularization_active_sink_count_(0),
        total_optimization_iterations_(0),
        objective_scale_promotion_count_(0),
        warned_regularization_budget_exceeded_(false) {
    validateOptions();
    initializeDecomposition();
  }

  DualDecomposition(int npartition, MinCutGraph min_cut_graph,
                    DualDecompositionOptions options = {})
      : DualDecomposition(npartition, min_cut_graph.nnode, min_cut_graph.narc,
                          std::move(min_cut_graph.arcs),
                          std::move(min_cut_graph.arc_capacities),
                          std::move(min_cut_graph.terminal_capacities),
                          options) {}

  long getTotalSolveLoopTime() const { return solve_loop_time_; }
  long getScale() const { return scale_; }
  double getBestLowerBound() const { return max_lower_bound_; }
  long getBestLowerBoundRaw() const { return max_lower_bound_raw_; }
  double getBestCertifiedLowerBound() const { return max_lower_bound_; }
  long getBestCertifiedLowerBoundRaw() const { return max_lower_bound_raw_; }
  double getBestRegularizedObjective() const {
    return max_regularized_objective_raw_ == std::numeric_limits<long>::min()
               ? -std::numeric_limits<double>::infinity()
               : double(max_regularized_objective_raw_) / scale_;
  }
  long getBestRegularizedObjectiveRaw() const {
    return max_regularized_objective_raw_;
  }
  long getBestUpperBoundRaw() const { return best_upper_bound_; }
  double getBestUpperBound() const {
    return best_upper_bound_ == std::numeric_limits<long>::max()
               ? std::numeric_limits<double>::infinity()
               : double(best_upper_bound_) / scale_;
  }
  long getCurrentUpperBoundRaw() const { return current_upper_bound_; }
  long getLastDisagreementCount() const { return last_disagreement_count_; }
  double getLastDisagreementNormSq() const {
    return last_disagreement_norm_sq_;
  }
  long getLastRegularizationBudget() const {
    return last_regularization_budget_;
  }
  long getLastRegularizationContribution() const {
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
  long getObjectiveScalePromotionCount() const {
    return objective_scale_promotion_count_;
  }
  const std::vector<PartitionPackage> &getPartitionPackages() const {
    if (!options_.emit_partition_packages) {
      throw std::runtime_error(
          "partition package export is disabled for this DualDecomposition");
    }
    return partition_packages_;
  }

  int regularizationStrengthForStepSize(long step_size) const {
    if (options_.regularization_scheme !=
        DualDecompositionRegularizationScheme::SCALED_EPSILON) {
      return 0;
    }
    return step_size <= 10 ? static_cast<int>(step_size) : 0;
  }

  void runPrimalSolutionDecodingStep(bool do_narrow_band_decode = false) {
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
      printf("recalc primal: %ld \n", primal_solver_->getMinCutValue());
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
    const int scaling_factor = 10;
    long step_size = options_.initial_step_size;
    scale_ = options_.objective_scale;
    total_optimization_iterations_ = 0;
    int iscale = 0;
    while (iscale < options_.num_optimization_scales && step_size >= 1) {
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
      if (status == REGULARIZATION_BUDGET_EXCEEDED &&
          tryPromoteObjectiveScale(/*factor=*/10, &step_size)) {
        iscale = 0;
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
      step_size /= 10;
      ++iscale;
      //auto rescale_problem_time =
      //    time_lambda([&] { scaleProblem<scaling_factor>(); });
      //printf("rescale problem time: %lums\n", rescale_problem_time.count());
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
    OptimizationStatus opt_status = ITERATION_COUNT_EXCEEDED;
    const bool report_progress = dualdecomp_progress_enabled();
    const auto scale_start = std::chrono::steady_clock::now();
    const int num_stats_in_group = 10;
    TwoGroupScalarStatisticsTracker<long> lower_bound_group_stats(
        num_stats_in_group);
    CycleCountingList dual_cycle_list;
    long max_lower_bound = std::numeric_limits<long>::min();
    int last_improvement_iter = 0;
    for (auto &solver_uptr : solvers_) {
      int regularization_str = regularizationStrengthForStepSize(step_size);
      solver_uptr->setRegularizationStrength(regularization_str);
    }
    for (int i = 0; i < nstep; ++i) {
      ++total_optimization_iterations_;

      std::vector<long> lower_bound_terms(solvers_.size(), 0);
      std::vector<long> regularization_budget_terms(solvers_.size(), 0);
      std::vector<long> regularization_contribution_terms(solvers_.size(), 0);
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
      long original_objective =
          std::accumulate(lower_bound_terms.begin(), lower_bound_terms.end(),
                          static_cast<long>(0));
      last_regularization_budget_ =
          std::accumulate(regularization_budget_terms.begin(),
                          regularization_budget_terms.end(),
                          static_cast<long>(0));
      last_regularization_contribution_ =
          std::accumulate(regularization_contribution_terms.begin(),
                          regularization_contribution_terms.end(),
                          static_cast<long>(0));
      last_regularization_anchor_sink_count_ =
          std::accumulate(regularization_anchor_count_terms.begin(),
                          regularization_anchor_count_terms.end(),
                          static_cast<long>(0));
      last_regularization_active_sink_count_ =
          std::accumulate(regularization_active_count_terms.begin(),
                          regularization_active_count_terms.end(),
                          static_cast<long>(0));
      const long regularized_objective =
          regularizedObjectiveRaw(original_objective,
                                  last_regularization_contribution_);
      const long lower_bound = certifiedOriginalLowerBoundRaw(
          original_objective, last_regularization_contribution_,
          last_regularization_budget_);
      warnIfRegularizationBudgetExceeded(last_regularization_budget_,
                                         regularizationStrengthForStepSize(
                                             step_size));
      if (isRegularizationBudgetExceeded(last_regularization_budget_,
                                         regularizationStrengthForStepSize(
                                             step_size))) {
        if (report_progress) {
          std::fprintf(stderr,
                       "mcpd3_progress stage=dd_solve_stop "
                       "reason=regularization_budget_exceeded iter=%d "
                       "budget=%ld limit=%ld step_size=%ld scale=%ld\n",
                       i, last_regularization_budget_,
                       regularizationBudgetLimit(), step_size, scale_);
          std::fflush(stderr);
        }
        opt_status = REGULARIZATION_BUDGET_EXCEEDED;
        break;
      }
      if (options_.track_primal_upper_bound) {
        current_upper_bound_ = updatePrimalUpperBound();
      }

      LagrangeUpdateStats update_stats;
      auto lagrange_update_time = time_lambda([&] {
        update_stats =
            runLagrangeMultipliersUpdateStep(step_size, use_momentum,
                                             lower_bound);
        disagreeing_global_indices_ =
            std::move(update_stats.disagreeing_global_indices);
          });
      last_disagreement_count_ = update_stats.disagreement_count;
      last_disagreement_norm_sq_ = update_stats.disagreement_norm_sq;

      const int regularization_strength =
          regularizationStrengthForStepSize(step_size);
      if (report_progress) {
        const long best_lower_bound =
            std::max(max_lower_bound, lower_bound);
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
            double(lower_bound) / scale_, double(best_lower_bound) / scale_,
            double(lower_bound) / scale_, double(best_lower_bound) / scale_,
            double(regularized_objective) / scale_,
            double(std::max(max_regularized_objective_raw_,
                            regularized_objective)) /
                scale_,
            current_upper_bound_ == std::numeric_limits<long>::max()
                ? std::numeric_limits<double>::infinity()
                : double(current_upper_bound_) / scale_,
            current_upper_bound_ == std::numeric_limits<long>::max()
                ? std::numeric_limits<double>::infinity()
                : double(current_upper_bound_ - lower_bound) / scale_,
            static_cast<long>(disagreeing_global_indices_.size()),
            update_stats.disagreement_norm_sq, step_size,
            update_stats.effective_step_size, regularization_strength,
            double(last_regularization_budget_) / scale_,
            double(last_regularization_contribution_) / scale_,
            last_regularization_anchor_sink_count_,
            last_regularization_active_sink_count_,
            i - last_improvement_iter, solve_loop_time.count(),
            lagrange_update_time.count(), elapsed, eta);
        std::fflush(stderr);
      }
      if (options_.verbose) {
        printf("iter : %6d lower_bound : %8.6lf best_lower_bound : %8.6lf upper_bound : %8.6lf gap : %8.6lf num_disagreeing : %6ld disagreement_norm_sq : %8.1lf step_size : %8ld regularization_strength : %6d regularization_budget : %8.6lf regularization_contribution : %8.6lf regularization_anchor_sink_count : %6ld regularization_active_sink_count : %6ld iters_since_improvement : %6d solve_loop_time: %8ldms lagrange_update_time: %8ldms\n",
               i, double(lower_bound) / scale_,
               double(std::max(max_lower_bound, lower_bound)) / scale_,
               current_upper_bound_ == std::numeric_limits<long>::max()
                   ? std::numeric_limits<double>::infinity()
                   : double(current_upper_bound_) / scale_,
               current_upper_bound_ == std::numeric_limits<long>::max()
                   ? std::numeric_limits<double>::infinity()
                   : double(current_upper_bound_ - lower_bound) / scale_,
               disagreeing_global_indices_.size(),
               update_stats.disagreement_norm_sq, update_stats.effective_step_size,
               regularization_strength,
               double(last_regularization_budget_) / scale_,
               double(last_regularization_contribution_) / scale_,
               last_regularization_anchor_sink_count_,
               last_regularization_active_sink_count_,
               i - last_improvement_iter, solve_loop_time.count(),
               lagrange_update_time.count());
      }

      max_lower_bound_ =
          std::max<double>(max_lower_bound_, double(lower_bound) / scale_);
      max_lower_bound_raw_ = std::max(max_lower_bound_raw_, lower_bound);
      max_regularized_objective_raw_ =
          std::max(max_regularized_objective_raw_, regularized_objective);

      if ( lower_bound > max_lower_bound ) {
        max_lower_bound = lower_bound;
        if (options_.legacy_patience &&
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
      } else if (!options_.legacy_patience &&
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

      if (best_upper_bound_ != std::numeric_limits<long>::max() &&
          max_lower_bound >= best_upper_bound_) {
        if (report_progress) {
          if (regularization_strength == 0) {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=lower_bound_closed_upper lower=%.6lf "
                         "upper=%.6lf\n",
                         double(max_lower_bound) / scale_,
                         double(best_upper_bound_) / scale_);
          } else {
            std::fprintf(stderr,
                         "mcpd3_progress stage=dd_solve_stop "
                         "reason=regularized_closed_upper lower=%.6lf "
                         "upper=%.6lf regularization_strength=%d\n",
                         double(max_lower_bound) / scale_,
                         double(best_upper_bound_) / scale_,
                         regularization_strength);
          }
          std::fflush(stderr);
        }
        if (options_.verbose) {
          if (regularization_strength == 0) {
            printf("breaking because lower bound closed primal upper bound: lower=%8.6lf upper=%8.6lf\n",
                   double(max_lower_bound) / scale_,
                   double(best_upper_bound_) / scale_);
          } else {
            printf("breaking because scaled epsilon regularization closed primal upper bound: lower=%8.6lf upper=%8.6lf regularization_strength=%d\n",
                   double(max_lower_bound) / scale_,
                   double(best_upper_bound_) / scale_,
                   regularization_strength);
          }
        }
        opt_status = OPTIMAL;
        break;
      }

      lower_bound_group_stats.addValue(lower_bound);
      if (options_.enable_group_stopping &&
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
                         i, double(lower_bound) / scale_);
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
                         double(lower_bound) / scale_);
            std::fflush(stderr);
          }
          if (options_.verbose) {
            printf("breaking because scaled epsilon regularized subproblems agree: regularization_strength=%d regularization_budget=%ld regularization_budget_limit=%ld\n",
                   regularization_strength, last_regularization_budget_,
                   regularizationBudgetLimit());
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
    }
    if (options_.verbose) {
      printf(" === MAX === lower_bound : %ld\n", max_lower_bound);
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
                   double(max_lower_bound) / scale_, elapsed);
      std::fflush(stderr);
    }
    return opt_status;
  }

  template <int scale> void scaleProblem() {
    scaleProblem(static_cast<long>(scale));
  }

  void scaleProblem(long scale) {
    if (scale <= 0) {
      throw std::runtime_error("problem scale factor must be positive");
    }
    scale_ = checkedScaleLong(scale_, scale);
    options_.objective_scale = checkedScaleLong(options_.objective_scale, scale);
    for (auto &cap : original_arc_capacities_) {
      cap = checkedScaleInt(cap, scale, options_.saturate_capacity_overflow);
    }
    for (auto &cap : original_terminal_capacities_) {
      cap = checkedScaleInt(cap, scale, options_.saturate_capacity_overflow);
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
        constraint.alpha = checkedScaleLong(constraint.alpha, scale);
        constraint.last_alpha = checkedScaleLong(constraint.last_alpha, scale);
      }
    }
    if (max_lower_bound_raw_ != std::numeric_limits<long>::min()) {
      max_lower_bound_raw_ = checkedScaleLong(max_lower_bound_raw_, scale);
    }
    if (max_regularized_objective_raw_ != std::numeric_limits<long>::min()) {
      max_regularized_objective_raw_ =
          checkedScaleLong(max_regularized_objective_raw_, scale);
    }
    if (best_upper_bound_ != std::numeric_limits<long>::max()) {
      best_upper_bound_ = checkedScaleLong(best_upper_bound_, scale);
    }
    if (current_upper_bound_ != std::numeric_limits<long>::max()) {
      current_upper_bound_ = checkedScaleLong(current_upper_bound_, scale);
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

  static int checkedScaleInt(int value, long scale,
                             bool saturate_capacity_overflow = false) {
    const long result = checkedScaleLong(value, scale);
    if (result > std::numeric_limits<int>::max() ||
        result < std::numeric_limits<int>::min()) {
      if (saturate_capacity_overflow) {
        return result < 0 ? std::numeric_limits<int>::min()
                          : std::numeric_limits<int>::max();
      }
      throw std::overflow_error("objective scale promotion exceeds int");
    }
    return static_cast<int>(result);
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
                           long lower_bound) const {
    std::hash<long> hasher{};
    long hash = hasher(lower_bound);
    for (const auto &global_index : disagreeing_global_indices) {
      hash ^= hasher(global_index);
    }
    return hash;
  }

  LagrangeUpdateStats runLagrangeMultipliersUpdateStep(long step_size,
                                                       bool use_momentum,
                                                       long lower_bound) {
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
    stats.effective_step_size =
        std::clamp(step_size, options_.min_step_size, options_.max_step_size);

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
            constraint.alpha += alpha_update;
          } else {
            constraint.alpha += stats.effective_step_size * diff;
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
    if (options_.initial_alpha_random_radius < 0) {
      throw std::runtime_error(
          "initial alpha random radius must be non-negative");
    }
    if (options_.regularization_budget_limit < 0) {
      throw std::runtime_error(
          "regularization budget limit must be non-negative");
    }
    if (options_.max_objective_scale_promotions < 0) {
      throw std::runtime_error(
          "max objective scale promotions must be non-negative");
    }
  }

  long regularizationBudgetLimit() const {
    return options_.regularization_budget_limit > 0
               ? options_.regularization_budget_limit
               : options_.objective_scale;
  }

  void warnIfRegularizationBudgetExceeded(long budget,
                                          int regularization_strength) {
    if (!isRegularizationBudgetExceeded(budget, regularization_strength) ||
        warned_regularization_budget_exceeded_) {
      return;
    }
    std::fprintf(stderr,
                 "warning: regularization budget %ld is not below limit %ld; "
                 "a regularized agreement may not certify optimality\n",
                 budget, regularizationBudgetLimit());
    std::fflush(stderr);
    warned_regularization_budget_exceeded_ = true;
  }

  bool isRegularizationBudgetExceeded(long budget,
                                      int regularization_strength) const {
    return regularization_strength > 0 && budget >= regularizationBudgetLimit();
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

  long computePrimalCutValue(const std::vector<bool> &labels) const {
    long cut_value = 0;
    for (int i = 0; i < narc_; ++i) {
      const int s = original_arcs_[2 * i + 0];
      const int t = original_arcs_[2 * i + 1];
      const int forward_capacity = original_arc_capacities_[2 * i + 0];
      const int backward_capacity = original_arc_capacities_[2 * i + 1];
      if (!labels[s] && labels[t]) {
        cut_value += forward_capacity;
      } else if (labels[s] && !labels[t]) {
        cut_value += backward_capacity;
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      const int terminal_capacity = original_terminal_capacities_[i];
      if (!labels[i] && terminal_capacity < 0) {
        cut_value += -terminal_capacity;
      } else if (labels[i] && terminal_capacity > 0) {
        cut_value += terminal_capacity;
      }
    }
    return cut_value;
  }

  long updatePrimalUpperBound() {
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
    const long upper_bound = computePrimalCutValue(decoded);
    if (upper_bound < best_upper_bound_) {
      best_upper_bound_ = upper_bound;
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
    std::vector<int> partitions_;
    partitions_ = configured_graph_partition(npartition_, narc_, nnode_, arcs_,
                                             &arc_capacities_);
    validateAndReportPartition(partitions_);
    dualdecomp_progress_report("dd_partition_done", 1, 1, init_start);
    auto mapping_start = std::chrono::steady_clock::now();
    int mapping_done = 0;
    for (auto &min_cut_sub_graph : min_cut_sub_graphs_) {
      min_cut_sub_graph.initializeMapping(nnode_);
      ++mapping_done;
      dualdecomp_progress_report("dd_initialize_mapping", mapping_done,
                                 npartition_, mapping_start);
    }
    /**
     * step 1: distribute all arcs into one and only one sub graph
     */
    auto arc_start = std::chrono::steady_clock::now();
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
      if (report_progress && (i + 1) % progress_interval == 0) {
        dualdecomp_progress_report("dd_distribute_arcs", i + 1, narc_,
                                   arc_start);
      }
    }
    dualdecomp_progress_report("dd_distribute_arcs", narc_, narc_, arc_start);
    arcs_.clear();
    arcs_.shrink_to_fit();
    arc_capacities_.clear();
    arc_capacities_.shrink_to_fit();
    /**
     * step 2: add source and sink capacities of nodes
     */
    auto terminal_start = std::chrono::steady_clock::now();
    for (int i = 0; i < nnode_; ++i) {
      if (terminal_capacities_[i] == 0) {
        if (report_progress && (i + 1) % progress_interval == 0) {
          dualdecomp_progress_report("dd_distribute_terminals", i + 1, nnode_,
                                     terminal_start);
        }
        continue;
      }
      min_cut_sub_graphs_[partitions_[i]].insertTerminal(
          i, terminal_capacities_[i]);
      if (report_progress && (i + 1) % progress_interval == 0) {
        dualdecomp_progress_report("dd_distribute_terminals", i + 1, nnode_,
                                   terminal_start);
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
        package.arcs = min_cut_sub_graph.graph.arcs;
        package.arc_capacities = min_cut_sub_graph.graph.arc_capacities;
        package.terminal_capacities =
            min_cut_sub_graph.graph.terminal_capacities;
        package.local_to_global = min_cut_sub_graph.local_to_global;
        package.constraint_endpoints.clear();
      }

      solvers_.emplace_back(std::make_unique<PrimalDualMinCutSolver>(
          std::move(min_cut_sub_graph.graph)));
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
    auto initial_alpha = [&]() -> long {
      if (!options_.randomize_initial_alphas ||
          options_.initial_alpha_random_radius == 0) {
        return 0;
      }
      return initial_alpha_distribution(initial_alpha_generator);
    };
    for (auto &[global_index, partitions] : constrained_nodes) {
      partitions.insert(
          partitions_[global_index]); // list each constrained node in its
                                      // original partition
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
          const long alpha = initial_alpha();
          constraint_arcs.emplace_back(
              /*alpha=*/alpha,
              /*last_alpha=*/alpha,
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
    partitions_.clear();
    partitions_.shrink_to_fit();
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
  std::vector<int> arc_capacities_;
  std::vector<int> terminal_capacities_;
  std::vector<int> original_arcs_;
  std::vector<int> original_arc_capacities_;
  std::vector<int> original_terminal_capacities_;

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
                   int forward_capacity, int backward_capacity) {
      int s = getOrInsertNode(global_source_index);
      int t = getOrInsertNode(global_target_index);
      graph.arc_capacities.push_back(forward_capacity);
      graph.arc_capacities.push_back(backward_capacity);
      graph.arcs.push_back(s);
      graph.arcs.push_back(t);
      graph.narc++;
    }

    void insertTerminal(int global_index, int terminal_capacity) {
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
  double max_lower_bound_;
  long max_lower_bound_raw_;
  long max_regularized_objective_raw_;
  long best_upper_bound_;
  long current_upper_bound_;
  long last_disagreement_count_;
  double last_disagreement_norm_sq_;
  long last_regularization_budget_;
  long last_regularization_contribution_;
  long last_regularization_anchor_sink_count_;
  long last_regularization_active_sink_count_;
  long total_optimization_iterations_;
  long objective_scale_promotion_count_;
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
