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
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <decomp/partition_worker.h>

namespace mcpd3 {

struct PartitionWorkerCoordinatorOptions {
  long min_step_size = 1;
  long max_step_size = 10000;
  bool use_momentum = true;
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
