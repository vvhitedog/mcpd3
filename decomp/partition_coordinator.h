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
#include <cstdio>
#include <functional>
#include <future>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <decomp/lower_bound_certificate.h>
#include <decomp/partition_worker.h>

namespace mcpd3 {

enum class PartitionWorkerOptimizationStatus {
  OPTIMAL,
  NO_FURTHER_PROGRESS,
  ITERATION_COUNT_EXCEEDED,
  REGULARIZATION_BUDGET_EXCEEDED
};

enum class PartitionWorkerStopReason {
  NONE,
  NO_DISAGREEMENT,
  REGULARIZED_NO_DISAGREEMENT,
  ITERATION_COUNT_EXCEEDED,
  NO_LOWER_BOUND_IMPROVEMENT,
  LEGACY_PATIENCE,
  GROUP_STOPPING,
  REGULARIZATION_BUDGET_EXCEEDED
};

enum class PartitionWorkerRegularizationScheme {
  SCALED_EPSILON,
  NONE
};

struct PartitionWorkerProgressRecord;

struct PartitionWorkerCoordinatorOptions {
  int num_optimization_scales = 5;
  int max_iteration_count = 10000;
  long initial_step_size = 10000;
  int patience = 10;
  bool legacy_patience = false;
  long min_step_size = 1;
  long max_step_size = 10000;
  long objective_scale = 1;
  bool use_momentum = true;
  bool enable_group_stopping = true;
  PartitionWorkerRegularizationScheme regularization_scheme =
      PartitionWorkerRegularizationScheme::SCALED_EPSILON;
  long regularization_budget_limit = 0;
  bool promote_objective_scale_on_overbudget = true;
  int max_objective_scale_promotions = 4;
  bool randomize_initial_alphas = false;
  long initial_alpha_random_radius = 0;
  unsigned int initial_alpha_random_seed = 0;
  int progress_report_interval = 0;
  std::function<void(const PartitionWorkerProgressRecord &)> progress_callback;
};

struct PartitionWorkerRoundStats {
  long round_id = 0;
  long lower_bound = 0;
  long regularized_objective = 0;
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
  long regularized_objective = 0;
  long best_regularized_objective = 0;
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
  long best_regularized_objective_raw = std::numeric_limits<long>::min();
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
  long best_regularized_objective_raw = std::numeric_limits<long>::min();
  double best_regularized_objective =
      -std::numeric_limits<double>::infinity();
  long total_iterations = 0;
  long final_disagreement_count = 0;
  double final_disagreement_norm_sq = 0;
  long final_regularization_budget = 0;
  long final_regularization_contribution = 0;
  long final_regularization_anchor_sink_count = 0;
  long final_regularization_active_sink_count = 0;
  long objective_scale_promotion_count = 0;
  std::vector<PartitionWorkerProgressRecord> progress_records;
  std::vector<PartitionWorkerScaleResult> scale_results;
};

class PartitionWorkerCoordinator {
public:
  PartitionWorkerCoordinator(
      std::vector<PartitionPackage> packages,
      std::vector<std::unique_ptr<PartitionWorker>> workers,
      PartitionWorkerCoordinatorOptions options = {})
      : workers_(std::move(workers)), options_(options),
        warned_regularization_budget_exceeded_(false) {
    validateOptions();
    if (workers_.empty()) {
      throw std::runtime_error("at least one partition worker is required");
    }
    packages_.reserve(packages.size());
    const auto package_to_worker_index = assignPackagesToWorkers(packages);
    std::vector<bool> worker_has_partition(workers_.size(), false);
    for (size_t i = 0; i < packages.size(); ++i) {
      const auto &package = packages[i];
      const int partition_id = package.partition_id;
      if (partition_to_worker_index_.find(partition_id) !=
          partition_to_worker_index_.end()) {
        throw std::runtime_error("duplicate partition id " +
                                 std::to_string(partition_id));
      }
      const size_t worker_index = package_to_worker_index[i];
      partition_to_worker_index_.emplace(partition_id, worker_index);
      partition_to_package_index_.emplace(partition_id, packages_.size());
      if (!worker_has_partition[worker_index]) {
        worker_has_partition[worker_index] = true;
        active_worker_indices_.push_back(worker_index);
      }
      workers_[worker_index]->loadPartition(package);

      PartitionPackage coordinator_package;
      coordinator_package.partition_id = partition_id;
      coordinator_package.constraint_endpoints = package.constraint_endpoints;
      packages_.push_back(std::move(coordinator_package));
    }
    buildConstraints();
    dropCoordinatorPackagePayloads();
    randomizeInitialAlphas();
  }

  PartitionWorkerRoundStats runRound(long round_id, long scale, long step_size,
                                     int regularization_strength) {
    std::vector<std::vector<AlphaUpdate>> alpha_updates(packages_.size());
    std::vector<size_t> synced_constraint_indices;
    for (size_t constraint_index = 0; constraint_index < constraints_.size();
         ++constraint_index) {
      const auto &constraint = constraints_[constraint_index];
      if (!constraint.needs_sync) {
        continue;
      }
      AlphaUpdate update{constraint.constraint_id, constraint.alpha,
                         constraint.last_alpha, constraint.alpha_momentum};
      alpha_updates[packageIndexForPartition(
          constraint.source.partition_id)].push_back(update);
      alpha_updates[packageIndexForPartition(
          constraint.target.partition_id)].push_back(update);
      synced_constraint_indices.push_back(constraint_index);
    }

    std::vector<PartitionSolveResult> results(packages_.size());
    std::vector<std::vector<size_t>> packages_by_worker(workers_.size());
    for (size_t package_index = 0; package_index < packages_.size();
         ++package_index) {
      const int partition_id = packages_[package_index].partition_id;
      const size_t worker_index = workerIndexForPartition(partition_id);
      packages_by_worker[worker_index].push_back(package_index);
    }

    std::vector<std::future<void>> futures;
    futures.reserve(active_worker_indices_.size());
    for (const auto worker_index : active_worker_indices_) {
      futures.push_back(std::async(
          std::launch::async,
          [&, worker_index] {
            std::vector<PartitionSolveRequest> requests;
            requests.reserve(packages_by_worker[worker_index].size());
            for (const auto package_index : packages_by_worker[worker_index]) {
              const int partition_id = packages_[package_index].partition_id;
              PartitionSolveRequest request;
              request.round_id = round_id;
              request.partition_id = partition_id;
              request.scale = scale;
              request.regularization_strength = regularization_strength;
              request.alpha_updates = std::move(alpha_updates[package_index]);
              requests.push_back(std::move(request));
            }
            const auto worker_results =
                workers_[worker_index]->solveRoundBatch(requests);
            if (worker_results.size() != requests.size()) {
              throw std::runtime_error(
                  "worker batch result count does not match request count");
            }
            for (const auto &worker_result : worker_results) {
              const size_t package_index =
                  packageIndexForPartition(worker_result.partition_id);
              if (workerIndexForPartition(worker_result.partition_id) !=
                  worker_index) {
                throw std::runtime_error(
                    "worker returned a result for an unowned partition");
              }
              results[package_index] = worker_result;
            }
          }));
    }
    for (auto &future : futures) {
      future.get();
    }
    for (const auto constraint_index : synced_constraint_indices) {
      constraints_[constraint_index].needs_sync = false;
    }

    PartitionWorkerRoundStats stats;
    stats.round_id = round_id;
    stats.effective_step_size =
        std::clamp(step_size, options_.min_step_size, options_.max_step_size);
    gatherRoundTerms(results, &stats);
    warnIfRegularizationBudgetExceeded(stats.regularization_budget,
                                       regularization_strength);
    updateConstraintsFromLabels(
        results, &stats,
        !isRegularizationBudgetExceeded(stats.regularization_budget,
                                        regularization_strength));
    return stats;
  }

  PartitionWorkerCoordinatorSolveResult solve() {
    PartitionWorkerCoordinatorSolveResult result;
    result.scale = options_.objective_scale;

    long schedule_scale = options_.initial_step_size;
    long step_size = options_.initial_step_size;
    int scale_index = 0;
    while (scale_index < options_.num_optimization_scales && step_size >= 1) {
      auto scale_result =
          runOptimizationScale(schedule_scale, step_size, &result);
      result.scale_results.push_back(scale_result);
      result.status = scale_result.status;
      result.stop_reason = scale_result.stop_reason;
      if (scale_result.status ==
              PartitionWorkerOptimizationStatus::REGULARIZATION_BUDGET_EXCEEDED &&
          tryPromoteObjectiveScale(/*factor=*/10, &schedule_scale, &step_size,
                                   &result)) {
        scale_index = 0;
        continue;
      }
      if (scale_result.status == PartitionWorkerOptimizationStatus::OPTIMAL) {
        break;
      }
      schedule_scale /= 10;
      step_size /= 10;
      ++scale_index;
    }

    if (result.best_lower_bound_raw != std::numeric_limits<long>::min()) {
      result.best_lower_bound =
          static_cast<double>(result.best_lower_bound_raw) / result.scale;
    }
    if (result.best_regularized_objective_raw !=
        std::numeric_limits<long>::min()) {
      result.best_regularized_objective =
          static_cast<double>(result.best_regularized_objective_raw) /
          result.scale;
    }
    return result;
  }

private:
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
    if (options_.progress_report_interval < 0) {
      throw std::runtime_error(
          "progress report interval must be non-negative");
    }
  }

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
    bool needs_sync = false;
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
    auto record_round =
        [&](int iteration, int regularization_strength,
            const PartitionWorkerRoundStats &round_stats) {
          result->final_disagreement_count = round_stats.disagreement_count;
          result->final_disagreement_norm_sq =
              round_stats.disagreement_norm_sq;
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
          const long best_regularized_objective =
              std::max(scale_result.best_regularized_objective_raw,
                       round_stats.regularized_objective);
          PartitionWorkerProgressRecord record{
              /*scale=*/scale,
              /*iteration=*/iteration,
              /*total_iteration=*/result->total_iterations,
              /*max_iteration=*/options_.max_iteration_count,
              /*lower_bound=*/round_stats.lower_bound,
              /*best_lower_bound=*/best_lower_bound,
              /*regularized_objective=*/round_stats.regularized_objective,
              /*best_regularized_objective=*/best_regularized_objective,
              /*disagreement_count=*/round_stats.disagreement_count,
              /*disagreement_norm_sq=*/round_stats.disagreement_norm_sq,
              /*step_size=*/step_size,
              /*effective_step_size=*/round_stats.effective_step_size,
              /*regularization_strength=*/regularization_strength,
              /*regularization_budget=*/round_stats.regularization_budget,
              /*regularization_contribution=*/
              round_stats.regularization_contribution,
              /*regularization_anchor_sink_count=*/
              round_stats.regularization_anchor_sink_count,
              /*regularization_active_sink_count=*/
              round_stats.regularization_active_sink_count,
              /*iterations_since_improvement=*/
              iteration - last_improvement_iter};
          result->progress_records.push_back(record);
          reportProgress(record);

          result->best_lower_bound_raw =
              std::max(result->best_lower_bound_raw, round_stats.lower_bound);
          scale_result.best_lower_bound_raw =
              std::max(scale_result.best_lower_bound_raw,
                       round_stats.lower_bound);
          result->best_regularized_objective_raw =
              std::max(result->best_regularized_objective_raw,
                       round_stats.regularized_objective);
          scale_result.best_regularized_objective_raw =
              std::max(scale_result.best_regularized_objective_raw,
                       round_stats.regularized_objective);
        };

    for (int i = 0; i < options_.max_iteration_count; ++i) {
      ++result->total_iterations;
      ++scale_result.iterations;

      const int regularization_strength =
          localRegularizationStrength(step_size);
      const auto round_stats =
          runRound(result->total_iterations, scale, step_size,
                   regularization_strength);

      if (isRegularizationBudgetExceeded(round_stats.regularization_budget,
                                         regularization_strength)) {
        scale_result.status =
            PartitionWorkerOptimizationStatus::REGULARIZATION_BUDGET_EXCEEDED;
        scale_result.stop_reason =
            PartitionWorkerStopReason::REGULARIZATION_BUDGET_EXCEEDED;
        scale_result.final_disagreement_count =
            round_stats.disagreement_count;
        scale_result.final_disagreement_norm_sq =
            round_stats.disagreement_norm_sq;
        return scale_result;
      }

      record_round(i, regularization_strength, round_stats);

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
          scale_result.status = PartitionWorkerOptimizationStatus::OPTIMAL;
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

  void dropCoordinatorPackagePayloads() {
    for (auto &package : packages_) {
      package.local_node_count = 0;
      package.arcs.clear();
      package.arc_capacities.clear();
      package.terminal_capacities.clear();
      package.local_to_global.clear();
      package.constraint_endpoints.clear();
      package.arcs.shrink_to_fit();
      package.arc_capacities.shrink_to_fit();
      package.terminal_capacities.shrink_to_fit();
      package.local_to_global.shrink_to_fit();
      package.constraint_endpoints.shrink_to_fit();
    }
  }

  void reportProgress(const PartitionWorkerProgressRecord &record) const {
    if (!options_.progress_callback ||
        options_.progress_report_interval <= 0) {
      return;
    }
    if (record.total_iteration % options_.progress_report_interval != 0) {
      return;
    }
    options_.progress_callback(record);
  }

  void randomizeInitialAlphas() {
    if (!options_.randomize_initial_alphas) {
      return;
    }
    if (options_.initial_alpha_random_radius < 0) {
      throw std::runtime_error(
          "initial alpha random radius must be non-negative");
    }
    if (options_.initial_alpha_random_radius == 0) {
      return;
    }

    std::mt19937 generator(options_.initial_alpha_random_seed);
    std::uniform_int_distribution<long> distribution(
        -options_.initial_alpha_random_radius,
        options_.initial_alpha_random_radius);
    for (auto &constraint : constraints_) {
      const long offset = distribution(generator);
      if (offset == 0) {
        continue;
      }
      constraint.alpha += offset;
      constraint.last_alpha += offset;
      constraint.needs_sync = true;
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

  static double estimatePartitionWork(const PartitionPackage &package) {
    const double node_count =
        static_cast<double>(std::max(package.local_node_count, 0));
    const double arc_count =
        static_cast<double>(package.arcs.size()) / 2.0;
    const double boundary_count =
        static_cast<double>(package.constraint_endpoints.size());
    return std::max(1.0, node_count + 2.0 * arc_count +
                             8.0 * boundary_count);
  }

  static double workerCapacity(
      const PartitionWorkerResourceEstimate &resources, double max_ram_gb) {
    const double cpu_count =
        static_cast<double>(std::max(resources.cpu_count, 1));
    double ram_factor = 1.0;
    if (max_ram_gb > 0.0) {
      const double ram_gb =
          static_cast<double>(std::max(resources.ram_gb, 0L));
      ram_factor = 0.75 + 0.25 * std::min(ram_gb / max_ram_gb, 1.0);
    }
    return std::max(1.0e-9, cpu_count * ram_factor);
  }

  std::vector<size_t> assignPackagesToWorkers(
      const std::vector<PartitionPackage> &packages) const {
    std::vector<PartitionWorkerResourceEstimate> resources;
    resources.reserve(workers_.size());
    double max_ram_gb = 0.0;
    for (const auto &worker : workers_) {
      resources.push_back(worker->resourceEstimate());
      max_ram_gb =
          std::max(max_ram_gb,
                   static_cast<double>(std::max(resources.back().ram_gb, 0L)));
    }

    std::vector<double> capacities;
    capacities.reserve(workers_.size());
    for (const auto &resource : resources) {
      capacities.push_back(workerCapacity(resource, max_ram_gb));
    }

    std::vector<double> package_work(packages.size(), 1.0);
    std::vector<size_t> package_order;
    package_order.reserve(packages.size());
    for (size_t i = 0; i < packages.size(); ++i) {
      package_work[i] = estimatePartitionWork(packages[i]);
      package_order.push_back(i);
    }
    std::sort(package_order.begin(), package_order.end(),
              [&](size_t lhs, size_t rhs) {
                if (package_work[lhs] != package_work[rhs]) {
                  return package_work[lhs] > package_work[rhs];
                }
                return lhs < rhs;
              });

    std::vector<double> worker_load(workers_.size(), 0.0);
    std::vector<size_t> package_to_worker(packages.size(), 0);
    constexpr double kTieTolerance = 1.0e-9;
    for (const auto package_index : package_order) {
      const double work = package_work[package_index];
      size_t best_worker = 0;
      double best_score = std::numeric_limits<double>::infinity();
      double best_raw_load = std::numeric_limits<double>::infinity();
      for (size_t worker_index = 0; worker_index < workers_.size();
           ++worker_index) {
        const double candidate_load = worker_load[worker_index] + work;
        const double candidate_score =
            candidate_load / capacities[worker_index];
        const bool better_score =
            candidate_score + kTieTolerance < best_score;
        const bool equal_score =
            std::fabs(candidate_score - best_score) <= kTieTolerance;
        const bool better_tie =
            equal_score &&
            (candidate_load + kTieTolerance < best_raw_load ||
             (std::fabs(candidate_load - best_raw_load) <= kTieTolerance &&
              worker_index < best_worker));
        if (better_score || better_tie) {
          best_worker = worker_index;
          best_score = candidate_score;
          best_raw_load = candidate_load;
        }
      }
      package_to_worker[package_index] = best_worker;
      worker_load[best_worker] += work;
    }
    return package_to_worker;
  }

  size_t packageIndexForPartition(int partition_id) const {
    auto find_iter = partition_to_package_index_.find(partition_id);
    if (find_iter == partition_to_package_index_.end()) {
      throw std::runtime_error("unknown package partition id " +
                               std::to_string(partition_id));
    }
    return find_iter->second;
  }

  void gatherRoundTerms(const std::vector<PartitionSolveResult> &results,
                        PartitionWorkerRoundStats *stats) const {
    long selected_lower_bound = 0;
    for (const auto &result : results) {
      selected_lower_bound =
          checkedAddObjectiveRaw(selected_lower_bound, result.lower_bound,
                                 "round lower bound overflow");
      stats->regularization_budget =
          checkedAddObjectiveRaw(stats->regularization_budget,
                                 result.regularization_budget,
                                 "round regularization budget overflow");
      stats->regularization_contribution =
          checkedAddObjectiveRaw(stats->regularization_contribution,
                                 result.regularization_contribution,
                                 "round regularization contribution overflow");
      stats->regularization_anchor_sink_count =
          checkedAddObjectiveRaw(stats->regularization_anchor_sink_count,
                                 result.regularization_anchor_sink_count,
                                 "round regularization anchor count overflow");
      stats->regularization_active_sink_count =
          checkedAddObjectiveRaw(stats->regularization_active_sink_count,
                                 result.regularization_active_sink_count,
                                 "round regularization active count overflow");
    }
    stats->regularized_objective =
        regularizedObjectiveRaw(selected_lower_bound,
                                stats->regularization_contribution);
    stats->lower_bound =
        certifiedOriginalLowerBoundRaw(selected_lower_bound,
                                       stats->regularization_contribution,
                                       stats->regularization_budget);
  }

  void updateConstraintsFromLabels(
      const std::vector<PartitionSolveResult> &results,
      PartitionWorkerRoundStats *stats, bool update_alpha = true) {
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

      if (!update_alpha) {
        continue;
      }
      const long old_alpha = constraint.alpha;
      const long old_last_alpha = constraint.last_alpha;
      const float old_alpha_momentum = constraint.alpha_momentum;
      constraint.last_alpha = constraint.alpha;
      if (diff == 0) {
        if (constraint.alpha != old_alpha ||
            constraint.last_alpha != old_last_alpha ||
            constraint.alpha_momentum != old_alpha_momentum) {
          constraint.needs_sync = true;
        }
        continue;
      }
      if (options_.use_momentum) {
        const double beta = .85;
        const int momentum_scale = 10;
        constraint.alpha_momentum =
            beta * constraint.alpha_momentum * beta + (1 - beta) * diff;
        const long alpha_update =
            stats->effective_step_size *
            static_cast<int>(momentum_scale * constraint.alpha_momentum);
        constraint.alpha += alpha_update;
      } else {
        constraint.alpha += stats->effective_step_size * diff;
      }
      if (constraint.alpha != old_alpha ||
          constraint.last_alpha != old_last_alpha ||
          constraint.alpha_momentum != old_alpha_momentum) {
        constraint.needs_sync = true;
      }
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

  static int checkedScaleInt(int value, long scale) {
    const long result = checkedScaleLong(value, scale);
    if (result > std::numeric_limits<int>::max() ||
        result < std::numeric_limits<int>::min()) {
      throw std::overflow_error("objective scale promotion exceeds int");
    }
    return static_cast<int>(result);
  }

  void scalePackage(PartitionPackage *package, long factor) {
    for (auto &capacity : package->arc_capacities) {
      capacity = checkedScaleInt(capacity, factor);
    }
    for (auto &capacity : package->terminal_capacities) {
      capacity = checkedScaleInt(capacity, factor);
    }
    for (auto &binding : package->constraint_endpoints) {
      binding.alpha = checkedScaleLong(binding.alpha, factor);
      binding.last_alpha = checkedScaleLong(binding.last_alpha, factor);
    }
  }

  void scaleObjectiveState(long factor,
                           PartitionWorkerCoordinatorSolveResult *result) {
    options_.objective_scale = checkedScaleLong(options_.objective_scale, factor);
    result->scale = options_.objective_scale;
    if (result->best_lower_bound_raw != std::numeric_limits<long>::min()) {
      result->best_lower_bound_raw =
          checkedScaleLong(result->best_lower_bound_raw, factor);
    }
    if (result->best_regularized_objective_raw !=
        std::numeric_limits<long>::min()) {
      result->best_regularized_objective_raw =
          checkedScaleLong(result->best_regularized_objective_raw, factor);
    }
    result->final_regularization_budget =
        checkedScaleLong(result->final_regularization_budget, factor);
    result->final_regularization_contribution =
        checkedScaleLong(result->final_regularization_contribution, factor);
    for (auto &package : packages_) {
      scalePackage(&package, factor);
    }
    for (auto &constraint : constraints_) {
      constraint.alpha = checkedScaleLong(constraint.alpha, factor);
      constraint.last_alpha = checkedScaleLong(constraint.last_alpha, factor);
    }
    for (const auto worker_index : active_worker_indices_) {
      workers_[worker_index]->scaleObjective(factor);
    }
    warned_regularization_budget_exceeded_ = false;
  }

  bool tryPromoteObjectiveScale(
      long factor, long *schedule_scale, long *step_size,
      PartitionWorkerCoordinatorSolveResult *result) {
    if (!options_.promote_objective_scale_on_overbudget ||
        options_.regularization_budget_limit > 0 ||
        result->objective_scale_promotion_count >=
            options_.max_objective_scale_promotions) {
      return false;
    }
    scaleObjectiveState(factor, result);
    ++result->objective_scale_promotion_count;
    *schedule_scale = options_.objective_scale;
    *step_size = options_.objective_scale;
    return true;
  }

  int localRegularizationStrength(long step_size) const {
    if (options_.regularization_scheme !=
        PartitionWorkerRegularizationScheme::SCALED_EPSILON) {
      return 0;
    }
    return step_size <= 10 ? static_cast<int>(step_size) : 0;
  }

  std::vector<PartitionPackage> packages_;
  std::vector<std::unique_ptr<PartitionWorker>> workers_;
  PartitionWorkerCoordinatorOptions options_;
  bool warned_regularization_budget_exceeded_;
  std::vector<size_t> active_worker_indices_;
  std::unordered_map<int, size_t> partition_to_worker_index_;
  std::unordered_map<int, size_t> partition_to_package_index_;
  std::vector<CoordinatorConstraint> constraints_;
  std::unordered_map<int, size_t> constraint_index_by_id_;
};

} // namespace mcpd3
