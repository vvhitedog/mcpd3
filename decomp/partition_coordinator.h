// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <decomp/partition_worker.h>

namespace mcpd3 {

enum class PartitionWorkerOptimizationStatus {
  OPTIMAL,
  NO_FURTHER_PROGRESS,
  ITERATION_COUNT_EXCEEDED
};

enum class PartitionWorkerStopReason {
  NONE,
  NO_DISAGREEMENT,
  REGULARIZED_NO_DISAGREEMENT,
  ITERATION_COUNT_EXCEEDED,
  NO_LOWER_BOUND_IMPROVEMENT,
  LEGACY_PATIENCE,
  GROUP_STOPPING
};

struct PartitionWorkerCoordinatorOptions {
  int num_optimization_scales = 5;
  int max_iteration_count = 10000;
  long initial_step_size = 10000;
  int patience = 10;
  bool legacy_patience = false;
  long min_step_size = 1;
  long max_step_size = 10000;
  bool use_momentum = true;
  bool enable_group_stopping = true;
};

struct PartitionWorkerRoundStats {
  long round_id = 0;
  long lower_bound = 0;
  long regularization_budget = 0;
  long regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
  long disagreement_count = 0;
  double disagreement_norm_sq = 0;
  long effective_step_size = 0;
  std::vector<int> disagreeing_global_indices;
};

struct PartitionWorkerProgressRecord {
  long scale = 1;
  int iteration = 0;
  long total_iteration = 0;
  int max_iteration = 0;
  long lower_bound = 0;
  long best_lower_bound = 0;
  long disagreement_count = 0;
  double disagreement_norm_sq = 0;
  long step_size = 0;
  long effective_step_size = 0;
  int regularization_strength = 0;
  long regularization_budget = 0;
  long regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
  int iterations_since_improvement = 0;
};

struct PartitionWorkerScaleResult {
  PartitionWorkerOptimizationStatus status =
      PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED;
  PartitionWorkerStopReason stop_reason =
      PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED;
  long scale = 1;
  long step_size = 0;
  int iterations = 0;
  long best_lower_bound_raw = std::numeric_limits<long>::min();
  long final_disagreement_count = 0;
  double final_disagreement_norm_sq = 0;
};

struct PartitionWorkerCoordinatorSolveResult {
  PartitionWorkerOptimizationStatus status =
      PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED;
  PartitionWorkerStopReason stop_reason =
      PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED;
  long scale = 1;
  long best_lower_bound_raw = std::numeric_limits<long>::min();
  double best_lower_bound = -std::numeric_limits<double>::infinity();
  long total_iterations = 0;
  long final_disagreement_count = 0;
  double final_disagreement_norm_sq = 0;
  long final_regularization_budget = 0;
  long final_regularization_contribution = 0;
  long final_regularization_anchor_sink_count = 0;
  long final_regularization_active_sink_count = 0;
  std::vector<PartitionWorkerProgressRecord> progress_records;
  std::vector<PartitionWorkerScaleResult> scale_results;
};

class PartitionWorkerCoordinator {
public:
  PartitionWorkerCoordinator(
      const std::vector<PartitionPackage> &packages,
      std::vector<std::unique_ptr<PartitionWorker>> workers,
      PartitionWorkerCoordinatorOptions options = {})
      : packages_(packages), workers_(std::move(workers)), options_(options) {
    if (packages_.size() != workers_.size()) {
      throw std::runtime_error(
          "partition package count must match worker count");
    }
    for (size_t i = 0; i < packages_.size(); ++i) {
      const int partition_id = packages_[i].partition_id;
      if (partition_to_worker_index_.find(partition_id) !=
          partition_to_worker_index_.end()) {
        throw std::runtime_error("duplicate partition id " +
                                 std::to_string(partition_id));
      }
      partition_to_worker_index_.emplace(partition_id, i);
      workers_[i]->loadPartition(packages_[i]);
    }
    buildConstraints();
  }

  PartitionWorkerRoundStats runRound(long round_id, long scale, long step_size,
                                     int regularization_strength) {
    std::vector<std::vector<AlphaUpdate>> alpha_updates(workers_.size());
    for (const auto &constraint : constraints_) {
      AlphaUpdate update{constraint.constraint_id, constraint.alpha,
                         constraint.last_alpha, constraint.alpha_momentum};
      alpha_updates[workerIndexForPartition(
          constraint.source.partition_id)].push_back(update);
      alpha_updates[workerIndexForPartition(
          constraint.target.partition_id)].push_back(update);
    }

    std::vector<PartitionSolveResult> results;
    results.reserve(workers_.size());
    for (size_t worker_index = 0; worker_index < workers_.size();
         ++worker_index) {
      PartitionSolveRequest request;
      request.round_id = round_id;
      request.scale = scale;
      request.regularization_strength = regularization_strength;
      request.alpha_updates = std::move(alpha_updates[worker_index]);
      results.push_back(workers_[worker_index]->solveRound(request));
    }

    PartitionWorkerRoundStats stats;
    stats.round_id = round_id;
    stats.effective_step_size =
        std::clamp(step_size, options_.min_step_size, options_.max_step_size);
    gatherRoundTerms(results, &stats);
    updateConstraintsFromLabels(results, &stats);
    return stats;
  }

  PartitionWorkerCoordinatorSolveResult solve() {
    PartitionWorkerCoordinatorSolveResult result;
    result.scale = options_.initial_step_size;

    long step_size = options_.initial_step_size;
    for (int scale_index = 0;
         scale_index < options_.num_optimization_scales; ++scale_index) {
      auto scale_result =
          runOptimizationScale(result.scale, step_size, &result);
      result.scale_results.push_back(scale_result);
      result.status = scale_result.status;
      result.stop_reason = scale_result.stop_reason;
      if (scale_result.status == PartitionWorkerOptimizationStatus::OPTIMAL) {
        break;
      }
      step_size /= 10;
    }

    if (result.best_lower_bound_raw != std::numeric_limits<long>::min()) {
      result.best_lower_bound =
          static_cast<double>(result.best_lower_bound_raw) / result.scale;
    }
    return result;
  }

private:
  struct ConstraintEndpoint {
    int partition_id = -1;
    int global_node_id = -1;
    int local_index = -1;
  };

  struct CoordinatorConstraint {
    int constraint_id = -1;
    ConstraintEndpoint source;
    ConstraintEndpoint target;
    long alpha = 0;
    long last_alpha = 0;
    float alpha_momentum = 0;
  };

  struct ConstraintAccumulator {
    int constraint_id = -1;
    bool has_source = false;
    bool has_target = false;
    ConstraintEndpoint source;
    ConstraintEndpoint target;
    long alpha = 0;
    long last_alpha = 0;
    float alpha_momentum = 0;
  };

  struct ConstraintLabels {
    bool has_source = false;
    bool has_target = false;
    int source_label = 0;
    int target_label = 0;
  };

  class TwoGroupMaxTracker {
  public:
    explicit TwoGroupMaxTracker(size_t group_size)
        : group_size_(group_size), ready_(false),
          first_group_max_(std::numeric_limits<long>::min()),
          second_group_max_(std::numeric_limits<long>::min()) {}

    void addValue(long value) {
      values_.push_back(value);
      if (values_.size() < 2 * group_size_) {
        ready_ = false;
        return;
      }

      first_group_max_ = *std::max_element(values_.begin(),
                                           values_.begin() + group_size_);
      second_group_max_ =
          *std::max_element(values_.begin() + group_size_, values_.end());
      values_.clear();
      ready_ = true;
    }

    bool areGroupsPopulated() const { return ready_; }

    std::pair<long, long> getMaximums() const {
      return {first_group_max_, second_group_max_};
    }

  private:
    size_t group_size_;
    bool ready_;
    long first_group_max_;
    long second_group_max_;
    std::vector<long> values_;
  };

  PartitionWorkerScaleResult runOptimizationScale(
      long scale, long step_size, PartitionWorkerCoordinatorSolveResult *result) {
    PartitionWorkerScaleResult scale_result;
    scale_result.scale = scale;
    scale_result.step_size = step_size;

    const int num_stats_in_group = 10;
    TwoGroupMaxTracker lower_bound_group_stats(num_stats_in_group);
    long scale_best_lower_bound = std::numeric_limits<long>::min();
    int last_improvement_iter = 0;

    for (int i = 0; i < options_.max_iteration_count; ++i) {
      ++result->total_iterations;
      ++scale_result.iterations;

      const int regularization_strength =
          step_size <= 10 ? static_cast<int>(step_size) : 0;
      const auto round_stats =
          runRound(result->total_iterations, scale, step_size,
                   regularization_strength);

      result->final_disagreement_count = round_stats.disagreement_count;
      result->final_disagreement_norm_sq = round_stats.disagreement_norm_sq;
      result->final_regularization_budget =
          round_stats.regularization_budget;
      result->final_regularization_contribution =
          round_stats.regularization_contribution;
      result->final_regularization_anchor_sink_count =
          round_stats.regularization_anchor_sink_count;
      result->final_regularization_active_sink_count =
          round_stats.regularization_active_sink_count;
      scale_result.final_disagreement_count =
          round_stats.disagreement_count;
      scale_result.final_disagreement_norm_sq =
          round_stats.disagreement_norm_sq;

      const long best_lower_bound =
          std::max(scale_best_lower_bound, round_stats.lower_bound);
      result->progress_records.push_back(
          PartitionWorkerProgressRecord{/*scale=*/scale,
                                        /*iteration=*/i,
                                        /*total_iteration=*/
                                        result->total_iterations,
                                        /*max_iteration=*/
                                        options_.max_iteration_count,
                                        /*lower_bound=*/
                                        round_stats.lower_bound,
                                        /*best_lower_bound=*/
                                        best_lower_bound,
                                        /*disagreement_count=*/
                                        round_stats.disagreement_count,
                                        /*disagreement_norm_sq=*/
                                        round_stats.disagreement_norm_sq,
                                        /*step_size=*/step_size,
                                        /*effective_step_size=*/
                                        round_stats.effective_step_size,
                                        /*regularization_strength=*/
                                        regularization_strength,
                                        /*regularization_budget=*/
                                        round_stats.regularization_budget,
                                        /*regularization_contribution=*/
                                        round_stats
                                            .regularization_contribution,
                                        /*regularization_anchor_sink_count=*/
                                        round_stats
                                            .regularization_anchor_sink_count,
                                        /*regularization_active_sink_count=*/
                                        round_stats
                                            .regularization_active_sink_count,
                                        /*iterations_since_improvement=*/
                                        i - last_improvement_iter});

      result->best_lower_bound_raw =
          std::max(result->best_lower_bound_raw, round_stats.lower_bound);
      scale_result.best_lower_bound_raw =
          std::max(scale_result.best_lower_bound_raw, round_stats.lower_bound);

      if (round_stats.lower_bound > scale_best_lower_bound) {
        scale_best_lower_bound = round_stats.lower_bound;
        if (options_.legacy_patience &&
            i - last_improvement_iter >= options_.patience) {
          scale_result.status =
              PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS;
          scale_result.stop_reason =
              PartitionWorkerStopReason::LEGACY_PATIENCE;
          return scale_result;
        }
        last_improvement_iter = i;
      } else if (!options_.legacy_patience &&
                 i - last_improvement_iter >= options_.patience) {
        scale_result.status =
            PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS;
        scale_result.stop_reason =
            PartitionWorkerStopReason::NO_LOWER_BOUND_IMPROVEMENT;
        return scale_result;
      }

      lower_bound_group_stats.addValue(round_stats.lower_bound);
      if (options_.enable_group_stopping &&
          lower_bound_group_stats.areGroupsPopulated()) {
        auto [first_group_max, second_group_max] =
            lower_bound_group_stats.getMaximums();
        if (second_group_max <= first_group_max) {
          scale_result.status =
              PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS;
          scale_result.stop_reason =
              PartitionWorkerStopReason::GROUP_STOPPING;
          return scale_result;
        }
      }

      if (round_stats.disagreeing_global_indices.empty()) {
        if (regularization_strength == 0) {
          scale_result.status = PartitionWorkerOptimizationStatus::OPTIMAL;
          scale_result.stop_reason =
              PartitionWorkerStopReason::NO_DISAGREEMENT;
        } else {
          scale_result.status =
              PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS;
          scale_result.stop_reason =
              PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT;
        }
        return scale_result;
      }
    }

    scale_result.status =
        PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED;
    scale_result.stop_reason =
        PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED;
    return scale_result;
  }

  void buildConstraints() {
    std::map<int, ConstraintAccumulator> accumulators;
    for (const auto &package : packages_) {
      for (const auto &binding : package.constraint_endpoints) {
        auto &accumulator = accumulators[binding.constraint_id];
        if (accumulator.constraint_id == -1) {
          accumulator.constraint_id = binding.constraint_id;
          accumulator.alpha = binding.alpha;
          accumulator.last_alpha = binding.last_alpha;
          accumulator.alpha_momentum = binding.alpha_momentum;
        }

        ConstraintEndpoint endpoint{package.partition_id,
                                    binding.global_node_id,
                                    binding.local_index};
        if (binding.is_source) {
          if (accumulator.has_source) {
            throw std::runtime_error("duplicate source endpoint for constraint " +
                                     std::to_string(binding.constraint_id));
          }
          accumulator.source = endpoint;
          accumulator.has_source = true;
        } else {
          if (accumulator.has_target) {
            throw std::runtime_error("duplicate target endpoint for constraint " +
                                     std::to_string(binding.constraint_id));
          }
          accumulator.target = endpoint;
          accumulator.has_target = true;
        }
      }
    }

    constraints_.clear();
    constraint_index_by_id_.clear();
    for (const auto &entry : accumulators) {
      const auto &accumulator = entry.second;
      if (!accumulator.has_source || !accumulator.has_target) {
        throw std::runtime_error("constraint " +
                                 std::to_string(accumulator.constraint_id) +
                                 " does not have source and target endpoints");
      }
      CoordinatorConstraint constraint;
      constraint.constraint_id = accumulator.constraint_id;
      constraint.source = accumulator.source;
      constraint.target = accumulator.target;
      constraint.alpha = accumulator.alpha;
      constraint.last_alpha = accumulator.last_alpha;
      constraint.alpha_momentum = accumulator.alpha_momentum;
      constraint_index_by_id_[constraint.constraint_id] = constraints_.size();
      constraints_.push_back(constraint);
    }
  }

  size_t workerIndexForPartition(int partition_id) const {
    auto find_iter = partition_to_worker_index_.find(partition_id);
    if (find_iter == partition_to_worker_index_.end()) {
      throw std::runtime_error("unknown worker partition id " +
                               std::to_string(partition_id));
    }
    return find_iter->second;
  }

  void gatherRoundTerms(const std::vector<PartitionSolveResult> &results,
                        PartitionWorkerRoundStats *stats) const {
    for (const auto &result : results) {
      stats->lower_bound += result.lower_bound;
      stats->regularization_budget += result.regularization_budget;
      stats->regularization_contribution +=
          result.regularization_contribution;
      stats->regularization_anchor_sink_count +=
          result.regularization_anchor_sink_count;
      stats->regularization_active_sink_count +=
          result.regularization_active_sink_count;
    }
  }

  void updateConstraintsFromLabels(
      const std::vector<PartitionSolveResult> &results,
      PartitionWorkerRoundStats *stats) {
    std::vector<ConstraintLabels> labels(constraints_.size());
    for (const auto &result : results) {
      for (const auto &label : result.constrained_labels) {
        auto constraint_iter = constraint_index_by_id_.find(label.constraint_id);
        if (constraint_iter == constraint_index_by_id_.end()) {
          throw std::runtime_error("unknown result constraint id " +
                                   std::to_string(label.constraint_id));
        }
        const size_t constraint_index = constraint_iter->second;
        const auto &constraint = constraints_[constraint_index];
        auto &constraint_labels = labels[constraint_index];
        if (result.partition_id == constraint.source.partition_id) {
          constraint_labels.source_label = label.label;
          constraint_labels.has_source = true;
        } else if (result.partition_id == constraint.target.partition_id) {
          constraint_labels.target_label = label.label;
          constraint_labels.has_target = true;
        } else {
          throw std::runtime_error("constraint label came from wrong partition");
        }
      }
    }

    for (size_t i = 0; i < constraints_.size(); ++i) {
      auto &constraint = constraints_[i];
      const auto &constraint_labels = labels[i];
      if (!constraint_labels.has_source || !constraint_labels.has_target) {
        throw std::runtime_error("missing labels for constraint " +
                                 std::to_string(constraint.constraint_id));
      }

      const int diff =
          constraint_labels.target_label - constraint_labels.source_label;
      if (diff != 0) {
        stats->disagreement_count += std::abs(diff);
        stats->disagreement_norm_sq += static_cast<double>(diff * diff);
        stats->disagreeing_global_indices.push_back(
            constraint.source.global_node_id);
      }

      constraint.last_alpha = constraint.alpha;
      if (diff == 0) {
        continue;
      }
      if (options_.use_momentum) {
        const double beta = .85;
        const int momentum_scale = 10;
        constraint.alpha_momentum =
            beta * constraint.alpha_momentum * beta + (1 - beta) * diff;
        constraint.alpha +=
            stats->effective_step_size *
            static_cast<int>(momentum_scale * constraint.alpha_momentum);
      } else {
        constraint.alpha += stats->effective_step_size * diff;
      }
    }
  }

  std::vector<PartitionPackage> packages_;
  std::vector<std::unique_ptr<PartitionWorker>> workers_;
  PartitionWorkerCoordinatorOptions options_;
  std::unordered_map<int, size_t> partition_to_worker_index_;
  std::vector<CoordinatorConstraint> constraints_;
  std::unordered_map<int, size_t> constraint_index_by_id_;
};

} // namespace mcpd3
