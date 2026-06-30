// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <future>
#include <limits>
#include <list>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <decomp/constraint.h>
#include <primaldual/mcpd3.h>

namespace mcpd3 {

struct ConstraintEndpointBinding {
  int constraint_id = -1;
  int global_node_id = -1;
  int local_index = -1;
  bool is_source = true;
  long alpha = 0;
  long last_alpha = 0;
  float alpha_momentum = 0;
};

struct AlphaUpdate {
  int constraint_id = -1;
  long alpha = 0;
  long last_alpha = 0;
  float alpha_momentum = 0;
};

struct ConstraintLabel {
  int constraint_id = -1;
  int global_node_id = -1;
  int local_index = -1;
  int label = 0;
};

struct PartitionPackage {
  int partition_id = -1;
  int local_node_count = 0;
  std::vector<int> arcs;
  std::vector<int> arc_capacities;
  std::vector<int> terminal_capacities;
  std::vector<int> local_to_global;
  std::vector<ConstraintEndpointBinding> constraint_endpoints;
};

struct PartitionSolveRequest {
  long round_id = 0;
  int partition_id = -1;
  long scale = 1;
  int regularization_strength = 0;
  std::vector<AlphaUpdate> alpha_updates;
};

struct PartitionSolveResult {
  long round_id = 0;
  int partition_id = -1;
  long lower_bound = 0;
  long regularization_budget = 0;
  long regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
  std::vector<ConstraintLabel> constrained_labels;
};

class PartitionWorker {
public:
  virtual ~PartitionWorker() = default;
  virtual void loadPartition(const PartitionPackage &package) = 0;
  virtual PartitionSolveResult solveRound(
      const PartitionSolveRequest &request) = 0;
  virtual std::vector<PartitionSolveResult> solveRoundBatch(
      const std::vector<PartitionSolveRequest> &requests) {
    std::vector<PartitionSolveResult> results;
    results.reserve(requests.size());
    for (const auto &request : requests) {
      results.push_back(solveRound(request));
    }
    return results;
  }
  virtual void scaleObjective(long factor) = 0;
};

class InProcessPartitionWorker final : public PartitionWorker {
public:
  void loadPartition(const PartitionPackage &package) override {
    validatePackage(package);
    if (partitions_.find(package.partition_id) != partitions_.end()) {
      throw std::runtime_error("partition id " +
                               std::to_string(package.partition_id) +
                               " is already loaded");
    }

    std::vector<int> arcs = package.arcs;
    auto &loaded = partitions_[package.partition_id];
    loaded.package = package;
    loaded.solver = std::make_unique<PrimalDualMinCutSolver>(
        package.local_node_count, static_cast<int>(package.arcs.size() / 2),
        std::move(arcs), package.arc_capacities, package.terminal_capacities);

    for (const auto &binding : package.constraint_endpoints) {
      addConstraintEndpoint(&loaded, binding);
    }
  }

  PartitionSolveResult solveRound(
      const PartitionSolveRequest &request) override {
    if (partitions_.empty()) {
      throw std::runtime_error("partition must be loaded before solveRound");
    }
    auto &loaded = loadedPartitionForRequest(request);
    return solveLoadedPartition(&loaded, request);
  }

  std::vector<PartitionSolveResult> solveRoundBatch(
      const std::vector<PartitionSolveRequest> &requests) override {
    if (partitions_.empty()) {
      throw std::runtime_error(
          "partition must be loaded before solveRoundBatch");
    }
    std::vector<LoadedPartition *> loaded_partitions;
    loaded_partitions.reserve(requests.size());
    std::unordered_set<int> seen_partition_ids;
    for (const auto &request : requests) {
      auto &loaded = loadedPartitionForRequest(request);
      if (!seen_partition_ids.insert(loaded.package.partition_id).second) {
        throw std::runtime_error(
            "batch solve requests must target distinct partitions");
      }
      loaded_partitions.push_back(&loaded);
    }

    std::vector<PartitionSolveResult> results(requests.size());
    std::vector<std::future<void>> futures;
    futures.reserve(requests.size());
    for (size_t i = 0; i < requests.size(); ++i) {
      futures.push_back(std::async(std::launch::async, [&, i] {
        results[i] = solveLoadedPartition(loaded_partitions[i], requests[i]);
      }));
    }
    for (auto &future : futures) {
      future.get();
    }
    return results;
  }

  void scaleObjective(long factor) override {
    if (partitions_.empty()) {
      throw std::runtime_error("partition must be loaded before scaleObjective");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    for (auto &[partition_id, loaded] : partitions_) {
      (void)partition_id;
      for (auto &capacity : loaded.package.arc_capacities) {
        capacity = checkedScaleInt(capacity, factor);
      }
      for (auto &capacity : loaded.package.terminal_capacities) {
        capacity = checkedScaleInt(capacity, factor);
      }
      for (auto &binding : loaded.package.constraint_endpoints) {
        binding.alpha = checkedScaleLong(binding.alpha, factor);
        binding.last_alpha = checkedScaleLong(binding.last_alpha, factor);
      }
      for (auto &constraint_arc : loaded.constraint_arcs) {
        constraint_arc.alpha = checkedScaleLong(constraint_arc.alpha, factor);
        constraint_arc.last_alpha =
            checkedScaleLong(constraint_arc.last_alpha, factor);
      }
      loaded.solver->scaleProblem(factor);
    }
  }

private:
  struct LoadedPartition {
    PartitionPackage package;
    std::unique_ptr<PrimalDualMinCutSolver> solver;
    std::list<DualDecompositionConstraintArc> constraint_arcs;
    std::unordered_map<int, DualDecompositionConstraintArcReference>
        constraint_arc_by_id;
  };

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

  static void validatePackage(const PartitionPackage &package) {
    if (package.partition_id < 0) {
      throw std::runtime_error("partition id must be non-negative");
    }
    if (package.local_node_count < 0) {
      throw std::runtime_error("local node count must be non-negative");
    }
    if (package.arcs.size() % 2 != 0) {
      throw std::runtime_error("partition arcs must contain endpoint pairs");
    }
    if (package.arc_capacities.size() != package.arcs.size()) {
      throw std::runtime_error("arc capacity count must match arc endpoints");
    }
    if (package.terminal_capacities.size() !=
        static_cast<size_t>(package.local_node_count)) {
      throw std::runtime_error(
          "terminal capacity count must match local node count");
    }
    if (!package.local_to_global.empty() &&
        package.local_to_global.size() !=
            static_cast<size_t>(package.local_node_count)) {
      throw std::runtime_error(
          "local_to_global count must match local node count");
    }
    for (const auto &local_index : package.arcs) {
      if (local_index < 0 || local_index >= package.local_node_count) {
        throw std::runtime_error("arc endpoint is outside local node range");
      }
    }
    for (const auto &binding : package.constraint_endpoints) {
      if (binding.constraint_id < 0) {
        throw std::runtime_error("constraint id must be non-negative");
      }
      if (binding.local_index < 0 ||
          binding.local_index >= package.local_node_count) {
        throw std::runtime_error(
            "constraint endpoint is outside local node range");
      }
    }
  }

  LoadedPartition &loadedPartitionForRequest(
      const PartitionSolveRequest &request) {
    if (request.partition_id < 0) {
      if (partitions_.size() != 1) {
        throw std::runtime_error(
            "partition id is required when multiple partitions are loaded");
      }
      return partitions_.begin()->second;
    }
    auto find_iter = partitions_.find(request.partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown solve request partition id " +
                               std::to_string(request.partition_id));
    }
    return find_iter->second;
  }

  void addConstraintEndpoint(LoadedPartition *loaded,
                             const ConstraintEndpointBinding &binding) {
    if (loaded->constraint_arc_by_id.find(binding.constraint_id) !=
        loaded->constraint_arc_by_id.end()) {
      throw std::runtime_error("duplicate constraint endpoint id " +
                               std::to_string(binding.constraint_id));
    }

    const int source_partition =
        binding.is_source ? loaded->package.partition_id : -1;
    const int target_partition =
        binding.is_source ? -1 : loaded->package.partition_id;
    const int source_local_index = binding.is_source ? binding.local_index : -1;
    const int target_local_index = binding.is_source ? -1 : binding.local_index;
    loaded->constraint_arcs.emplace_back(
        binding.alpha, binding.last_alpha, binding.alpha_momentum,
        source_partition, target_partition, source_local_index,
        target_local_index);
    auto arc_reference = --loaded->constraint_arcs.end();
    loaded->constraint_arc_by_id.emplace(binding.constraint_id, arc_reference);
    if (binding.is_source) {
      loaded->solver->addSourceDualDecompositionConstraint(arc_reference);
    } else {
      loaded->solver->addTargetDualDecompositionConstraint(arc_reference);
    }
  }

  void applyAlphaUpdate(LoadedPartition *loaded, const AlphaUpdate &update) {
    auto find_iter = loaded->constraint_arc_by_id.find(update.constraint_id);
    if (find_iter == loaded->constraint_arc_by_id.end()) {
      throw std::runtime_error("unknown alpha update constraint id " +
                               std::to_string(update.constraint_id));
    }
    auto arc_reference = find_iter->second;
    arc_reference->alpha = update.alpha;
    arc_reference->last_alpha = update.last_alpha;
    arc_reference->alpha_momentum = update.alpha_momentum;
  }

  PartitionSolveResult solveLoadedPartition(
      LoadedPartition *loaded, const PartitionSolveRequest &request) {
    for (const auto &update : request.alpha_updates) {
      applyAlphaUpdate(loaded, update);
    }

    loaded->solver->setRegularizationStrength(request.regularization_strength);
    loaded->solver->solve();

    PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = loaded->package.partition_id;
    result.lower_bound = loaded->solver->getMinCutValue();
    result.regularization_budget =
        loaded->solver->getLastRegularizationBudget();
    result.regularization_contribution =
        loaded->solver->getLastRegularizationContribution();
    result.regularization_anchor_sink_count =
        loaded->solver->getLastRegularizationAnchorSinkCount();
    result.regularization_active_sink_count =
        loaded->solver->getLastRegularizationActiveSinkCount();

    result.constrained_labels.reserve(
        loaded->package.constraint_endpoints.size());
    for (const auto &binding : loaded->package.constraint_endpoints) {
      result.constrained_labels.push_back(ConstraintLabel{
          binding.constraint_id, binding.global_node_id, binding.local_index,
          loaded->solver->getMinCutSolution(binding.local_index)});
    }
    return result;
  }

  std::unordered_map<int, LoadedPartition> partitions_;
};

} // namespace mcpd3
