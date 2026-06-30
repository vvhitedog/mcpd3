// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <limits>
#include <list>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
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
  virtual void scaleObjective(long factor) = 0;
};

class InProcessPartitionWorker final : public PartitionWorker {
public:
  void loadPartition(const PartitionPackage &package) override {
    validatePackage(package);

    package_ = package;
    solver_.reset();
    constraint_arcs_.clear();
    constraint_arc_by_id_.clear();

    std::vector<int> arcs = package.arcs;
    solver_ = std::make_unique<PrimalDualMinCutSolver>(
        package.local_node_count, static_cast<int>(package.arcs.size() / 2),
        std::move(arcs), package.arc_capacities, package.terminal_capacities);

    for (const auto &binding : package.constraint_endpoints) {
      addConstraintEndpoint(binding);
    }
  }

  PartitionSolveResult solveRound(
      const PartitionSolveRequest &request) override {
    if (!solver_) {
      throw std::runtime_error("partition must be loaded before solveRound");
    }

    for (const auto &update : request.alpha_updates) {
      applyAlphaUpdate(update);
    }

    solver_->setRegularizationStrength(request.regularization_strength);
    solver_->solve();

    PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = package_.partition_id;
    result.lower_bound = solver_->getMinCutValue();
    result.regularization_budget = solver_->getLastRegularizationBudget();
    result.regularization_contribution =
        solver_->getLastRegularizationContribution();
    result.regularization_anchor_sink_count =
        solver_->getLastRegularizationAnchorSinkCount();
    result.regularization_active_sink_count =
        solver_->getLastRegularizationActiveSinkCount();

    result.constrained_labels.reserve(package_.constraint_endpoints.size());
    for (const auto &binding : package_.constraint_endpoints) {
      result.constrained_labels.push_back(ConstraintLabel{
          binding.constraint_id, binding.global_node_id, binding.local_index,
          solver_->getMinCutSolution(binding.local_index)});
    }
    return result;
  }

  void scaleObjective(long factor) override {
    if (!solver_) {
      throw std::runtime_error("partition must be loaded before scaleObjective");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    for (auto &capacity : package_.arc_capacities) {
      capacity = checkedScaleInt(capacity, factor);
    }
    for (auto &capacity : package_.terminal_capacities) {
      capacity = checkedScaleInt(capacity, factor);
    }
    for (auto &binding : package_.constraint_endpoints) {
      binding.alpha = checkedScaleLong(binding.alpha, factor);
      binding.last_alpha = checkedScaleLong(binding.last_alpha, factor);
    }
    for (auto &constraint_arc : constraint_arcs_) {
      constraint_arc.alpha = checkedScaleLong(constraint_arc.alpha, factor);
      constraint_arc.last_alpha =
          checkedScaleLong(constraint_arc.last_alpha, factor);
    }
    solver_->scaleProblem(factor);
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

  void addConstraintEndpoint(const ConstraintEndpointBinding &binding) {
    if (constraint_arc_by_id_.find(binding.constraint_id) !=
        constraint_arc_by_id_.end()) {
      throw std::runtime_error("duplicate constraint endpoint id " +
                               std::to_string(binding.constraint_id));
    }

    const int source_partition =
        binding.is_source ? package_.partition_id : -1;
    const int target_partition =
        binding.is_source ? -1 : package_.partition_id;
    const int source_local_index = binding.is_source ? binding.local_index : -1;
    const int target_local_index = binding.is_source ? -1 : binding.local_index;
    constraint_arcs_.emplace_back(binding.alpha, binding.last_alpha,
                                  binding.alpha_momentum, source_partition,
                                  target_partition, source_local_index,
                                  target_local_index);
    auto arc_reference = --constraint_arcs_.end();
    constraint_arc_by_id_.emplace(binding.constraint_id, arc_reference);
    if (binding.is_source) {
      solver_->addSourceDualDecompositionConstraint(arc_reference);
    } else {
      solver_->addTargetDualDecompositionConstraint(arc_reference);
    }
  }

  void applyAlphaUpdate(const AlphaUpdate &update) {
    auto find_iter = constraint_arc_by_id_.find(update.constraint_id);
    if (find_iter == constraint_arc_by_id_.end()) {
      throw std::runtime_error("unknown alpha update constraint id " +
                               std::to_string(update.constraint_id));
    }
    auto arc_reference = find_iter->second;
    arc_reference->alpha = update.alpha;
    arc_reference->last_alpha = update.last_alpha;
    arc_reference->alpha_momentum = update.alpha_momentum;
  }

  PartitionPackage package_;
  std::unique_ptr<PrimalDualMinCutSolver> solver_;
  std::list<DualDecompositionConstraintArc> constraint_arcs_;
  std::unordered_map<int, DualDecompositionConstraintArcReference>
      constraint_arc_by_id_;
};

} // namespace mcpd3
