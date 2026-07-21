// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <future>
#include <limits>
#include <list>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <decomp/constraint.h>
#include <primaldual/mcpd3.h>

namespace mcpd3 {

inline bool inprocess_worker_sequential_batch_enabled() {
  const char *value = std::getenv("MCPD3_INPROCESS_WORKER_SEQUENTIAL_BATCH");
  return value != nullptr && value[0] != '\0' && std::string(value) != "0";
}

struct ConstraintEndpointBinding {
  int constraint_id = -1;
  int global_node_id = -1;
  int local_index = -1;
  bool is_source = true;
  Lagrange alpha = 0;
  Lagrange last_alpha = 0;
  float alpha_momentum = 0;
};

struct AlphaUpdate {
  int constraint_id = -1;
  Lagrange alpha = 0;
  Lagrange last_alpha = 0;
  float alpha_momentum = 0;
};

struct ConstraintLabel {
  int constraint_id = -1;
  int global_node_id = -1;
  int local_index = -1;
  int label = 0;
};

struct NodeLabel {
  int global_node_id = -1;
  int local_index = -1;
  int label = 0;
};

struct PartitionPackage {
  int partition_id = -1;
  int local_node_count = 0;
  SolverArray<int> arcs;
  SolverArray<Capacity> arc_capacities;
  SolverArray<Capacity> terminal_capacities;
  SolverArray<int> local_to_global;
  std::vector<ConstraintEndpointBinding> constraint_endpoints;
  long objective_multiplier = 1;
  CanonicalCutSelection canonical_cut_selection =
      CanonicalCutSelection::SOLVER_DEFAULT;
  bool force_full_mincut_recompute = false;
  SolverArray<int> reference_cut_labels;
  ReferenceCutSelection reference_cut_selection =
      ReferenceCutSelection::CLOSEST_EXACT;
  long reference_cut_check_interval = 1;

  std::size_t fileBackedBytes() const {
    return arcs.fileBackedBytes() + arc_capacities.fileBackedBytes() +
           terminal_capacities.fileBackedBytes() +
           local_to_global.fileBackedBytes() +
           reference_cut_labels.fileBackedBytes();
  }
};

struct PartitionSolveRequest {
  long round_id = 0;
  int partition_id = -1;
  long scale = 1;
  Capacity regularization_strength = 0;
  bool return_full_labels = false;
  std::vector<AlphaUpdate> alpha_updates;
};

struct PartitionSolveResult {
  long round_id = 0;
  int partition_id = -1;
  Objective lower_bound = 0;
  Objective regularization_budget = 0;
  Objective regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
  std::vector<ConstraintLabel> constrained_labels;
  std::vector<NodeLabel> full_labels;
};

struct PartitionCapacityUpdate {
  int partition_id = -1;
  SolverArray<Capacity> arc_capacities;
  SolverArray<Capacity> terminal_capacities;
  bool preserve_flow_state = true;
  Objective flow_scale_numerator = 1;
  Objective flow_scale_denominator = 1;
};

struct PartitionWorkerResourceEstimate {
  int cpu_count = 1;
  long ram_gb = 0;
};

inline Objective checkedScaleWorkerObjective(const Objective &value,
                                              long scale) {
  return checked_scale(value, scale, "objective scale promotion overflow");
}

inline Lagrange checkedScaleWorkerLagrange(const Lagrange &value,
                                           long scale) {
  return checked_scale(value, scale, "lagrange scale promotion overflow");
}

inline Capacity checkedScaleWorkerCapacity(
    const Capacity &value, long scale,
    bool saturate_capacity_overflow = false) {
  return checked_scale_capacity(value, scale, saturate_capacity_overflow);
}

inline void validatePartitionPackage(const PartitionPackage &package) {
  if (package.partition_id < 0) {
    throw std::runtime_error("partition id must be non-negative");
  }
  if (package.local_node_count < 0) {
    throw std::runtime_error("local node count must be non-negative");
  }
  if (package.objective_multiplier <= 0) {
    throw std::runtime_error("partition objective multiplier must be positive");
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
  if (!package.reference_cut_labels.empty() &&
      package.reference_cut_labels.size() !=
          static_cast<size_t>(package.local_node_count)) {
    throw std::runtime_error(
        "reference cut label count must match local node count");
  }
  if (package.reference_cut_check_interval <= 0) {
    throw std::runtime_error("reference cut check interval must be positive");
  }
}

class PartitionWorker {
public:
  virtual ~PartitionWorker() = default;
  virtual PartitionWorkerResourceEstimate resourceEstimate() const {
    return {};
  }
  virtual void loadPartition(const PartitionPackage &package) = 0;
  virtual void loadPartition(PartitionPackage &&package) {
    loadPartition(static_cast<const PartitionPackage &>(package));
  }
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
  virtual void scaleObjective(long factor,
                              bool saturate_capacity_overflow = false) = 0;
  virtual void scaleObjectivePartitions(
      const std::vector<int> &partition_ids, long factor,
      bool saturate_capacity_overflow = false) {
    if (partition_ids.empty()) {
      throw std::runtime_error(
          "objective scaling requires at least one partition id");
    }
    scaleObjective(factor, saturate_capacity_overflow);
  }
  virtual void replacePartitionCapacities(
      const PartitionCapacityUpdate &update) {
    (void)update;
    throw std::runtime_error(
        "partition worker does not support capacity replacement");
  }
  virtual void replacePartitionCapacitiesFor(
      int target_partition_id, const PartitionCapacityUpdate &update) {
    if (target_partition_id != update.partition_id) {
      throw std::runtime_error(
          "partition worker does not support remapped capacity replacement");
    }
    replacePartitionCapacities(update);
  }
  virtual std::size_t fullLabelCount(int partition_id) const {
    (void)partition_id;
    throw std::runtime_error(
        "partition worker does not support bounded label recovery");
  }
  virtual void copyFullLabels(int partition_id, std::size_t offset,
                              NodeLabel *destination, std::size_t count) {
    (void)partition_id;
    (void)offset;
    (void)destination;
    (void)count;
    throw std::runtime_error(
        "partition worker does not support bounded label recovery");
  }
};

class InProcessPartitionWorker final : public PartitionWorker {
public:
  InProcessPartitionWorker() = default;
  explicit InProcessPartitionWorker(SolverStorageOptions storage_options)
      : storage_options_(std::move(storage_options)) {}

  void loadPartition(const PartitionPackage &package) override {
    PartitionPackage copy = package;
    loadPartition(std::move(copy));
  }

  void loadPartition(PartitionPackage &&package) override {
    validatePackage(package);
    if (partitions_.find(package.partition_id) != partitions_.end()) {
      throw std::runtime_error("partition id " +
                               std::to_string(package.partition_id) +
                               " is already loaded");
    }

    const int partition_id = package.partition_id;
    const int local_node_count = package.local_node_count;
    const int arc_count = static_cast<int>(package.arcs.size() / 2);
    auto &loaded = partitions_[package.partition_id];
    loaded.partition_id = partition_id;
    loaded.local_node_count = local_node_count;
    loaded.local_to_global = std::move(package.local_to_global)
                                 .rehome(storage_options_, "local_to_global");
    loaded.constraint_endpoints = std::move(package.constraint_endpoints);
    loaded.solver = std::make_unique<PrimalDualMinCutSolver>(
        local_node_count, arc_count, std::move(package.arcs),
        std::move(package.arc_capacities),
        std::move(package.terminal_capacities), storage_options_);
    loaded.solver->setCanonicalCutSelection(package.canonical_cut_selection);
    loaded.solver->setForceFullMinCutRecompute(
        package.force_full_mincut_recompute);
    if (!package.reference_cut_labels.empty()) {
      loaded.solver->setReferenceCutLabelsStorage(
          std::move(package.reference_cut_labels));
      loaded.solver->setReferenceCutSelection(
          package.reference_cut_selection);
      loaded.solver->setReferenceCutCheckInterval(
          package.reference_cut_check_interval);
    }

    for (const auto &binding : loaded.constraint_endpoints) {
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
      if (!seen_partition_ids.insert(loaded.partition_id).second) {
        throw std::runtime_error(
            "batch solve requests must target distinct partitions");
      }
      loaded_partitions.push_back(&loaded);
    }

    if (inprocess_worker_sequential_batch_enabled()) {
      std::vector<PartitionSolveResult> results;
      results.reserve(requests.size());
      for (size_t i = 0; i < requests.size(); ++i) {
        results.push_back(
            solveLoadedPartition(loaded_partitions[i], requests[i]));
      }
      return results;
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

  void scaleObjective(long factor,
                      bool saturate_capacity_overflow = false) override {
    if (partitions_.empty()) {
      throw std::runtime_error("partition must be loaded before scaleObjective");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    for (auto &[partition_id, loaded] : partitions_) {
      (void)partition_id;
      for (auto &constraint_arc : loaded.constraint_arcs) {
        constraint_arc.alpha =
            checkedScaleWorkerLagrange(constraint_arc.alpha, factor);
        constraint_arc.last_alpha =
            checkedScaleWorkerLagrange(constraint_arc.last_alpha, factor);
      }
      loaded.solver->scaleProblem(factor, saturate_capacity_overflow);
    }
  }

  void scaleObjectivePartitions(
      const std::vector<int> &partition_ids, long factor,
      bool saturate_capacity_overflow = false) override {
    if (partition_ids.empty()) {
      throw std::runtime_error(
          "objective scaling requires at least one partition id");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    std::unordered_set<int> seen;
    std::vector<LoadedPartition *> loaded_partitions;
    loaded_partitions.reserve(partition_ids.size());
    for (const int partition_id : partition_ids) {
      if (!seen.insert(partition_id).second) {
        throw std::runtime_error("duplicate objective scale partition id " +
                                 std::to_string(partition_id));
      }
      loaded_partitions.push_back(&loadedPartitionById(partition_id));
    }
    for (auto *loaded : loaded_partitions) {
      for (auto &constraint_arc : loaded->constraint_arcs) {
        constraint_arc.alpha =
            checkedScaleWorkerLagrange(constraint_arc.alpha, factor);
        constraint_arc.last_alpha =
            checkedScaleWorkerLagrange(constraint_arc.last_alpha, factor);
      }
      loaded->solver->scaleProblem(factor, saturate_capacity_overflow);
    }
  }

  void replacePartitionCapacities(
      const PartitionCapacityUpdate &update) override {
    replacePartitionCapacitiesFor(update.partition_id, update);
  }

  void replacePartitionCapacitiesFor(
      int target_partition_id,
      const PartitionCapacityUpdate &update) override {
    auto &loaded = loadedPartitionById(target_partition_id);
    loaded.solver->replaceProblemCapacities(
        update.arc_capacities, update.terminal_capacities,
        update.preserve_flow_state, update.flow_scale_numerator,
        update.flow_scale_denominator);
  }

  std::size_t fullLabelCount(int partition_id) const override {
    return static_cast<std::size_t>(
        loadedPartitionById(partition_id).local_node_count);
  }

  void copyFullLabels(int partition_id, std::size_t offset,
                      NodeLabel *destination, std::size_t count) override {
    const auto &loaded = loadedPartitionById(partition_id);
    const auto total = static_cast<std::size_t>(loaded.local_node_count);
    if (offset > total || count > total - offset) {
      throw std::runtime_error("full label range is outside partition");
    }
    if (count != 0 && destination == nullptr) {
      throw std::runtime_error("full label destination must not be null");
    }
    for (std::size_t index = 0; index < count; ++index) {
      const auto local_index = offset + index;
      const int global_node_id =
          loaded.local_to_global.empty()
              ? static_cast<int>(local_index)
              : loaded.local_to_global[local_index];
      destination[index] = NodeLabel{
          global_node_id, static_cast<int>(local_index),
          loaded.solver->getMinCutSolution(static_cast<int>(local_index))};
    }
  }

private:
  struct LoadedPartition {
    int partition_id = -1;
    int local_node_count = 0;
    SolverArray<int> local_to_global;
    std::vector<ConstraintEndpointBinding> constraint_endpoints;
    std::unique_ptr<PrimalDualMinCutSolver> solver;
    std::list<DualDecompositionConstraintArc> constraint_arcs;
    std::unordered_map<int, DualDecompositionConstraintArcReference>
        constraint_arc_by_id;
  };

public:
  std::vector<int> minCutSolution(int partition_id) const {
    const auto &loaded = loadedPartitionById(partition_id);
    std::vector<int> solution;
    solution.reserve(static_cast<size_t>(loaded.local_node_count));
    for (int i = 0; i < loaded.local_node_count; ++i) {
      solution.push_back(loaded.solver->getMinCutSolution(i));
    }
    return solution;
  }

  void restoreMinCutSolution(int partition_id,
                             const std::vector<int> &solution) {
    auto &loaded = loadedPartitionById(partition_id);
    if (solution.size() != static_cast<size_t>(loaded.local_node_count)) {
      throw std::runtime_error("restored min-cut solution has wrong size");
    }
    loaded.solver->setMinCutSolution(solution);
  }

  PrimalDualMinCutSolver::WarmState warmState(int partition_id) const {
    const auto &loaded = loadedPartitionById(partition_id);
    return loaded.solver->captureWarmState();
  }

  PrimalDualMinCutSolver::StorageDiagnostics
  storageDiagnostics(int partition_id) const {
    return loadedPartitionById(partition_id).solver->getStorageDiagnostics();
  }

  void restoreWarmState(int partition_id,
                        const PrimalDualMinCutSolver::WarmState &state) {
    auto &loaded = loadedPartitionById(partition_id);
    loaded.solver->restoreWarmState(state);
  }

private:
  static void validatePackage(const PartitionPackage &package) {
    validatePartitionPackage(package);
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

  LoadedPartition &loadedPartitionById(int partition_id) {
    if (partition_id < 0) {
      if (partitions_.size() != 1) {
        throw std::runtime_error(
            "partition id is required when multiple partitions are loaded");
      }
      return partitions_.begin()->second;
    }
    auto find_iter = partitions_.find(partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown partition id " +
                               std::to_string(partition_id));
    }
    return find_iter->second;
  }

  const LoadedPartition &loadedPartitionById(int partition_id) const {
    if (partition_id < 0) {
      if (partitions_.size() != 1) {
        throw std::runtime_error(
            "partition id is required when multiple partitions are loaded");
      }
      return partitions_.begin()->second;
    }
    auto find_iter = partitions_.find(partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown partition id " +
                               std::to_string(partition_id));
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

    const int source_partition = binding.is_source ? loaded->partition_id : -1;
    const int target_partition = binding.is_source ? -1 : loaded->partition_id;
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
    arc_reference->last_alpha = arc_reference->alpha;
    arc_reference->alpha = update.alpha;
  }

  void markAlphaStateSolved(LoadedPartition *loaded) {
    for (auto &constraint_arc : loaded->constraint_arcs) {
      constraint_arc.last_alpha = constraint_arc.alpha;
    }
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
    result.partition_id = loaded->partition_id;
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
        loaded->constraint_endpoints.size());
    for (const auto &binding : loaded->constraint_endpoints) {
      result.constrained_labels.push_back(ConstraintLabel{
          binding.constraint_id, binding.global_node_id, binding.local_index,
          loaded->solver->getMinCutSolution(binding.local_index)});
    }
    if (request.return_full_labels) {
      result.full_labels.reserve(static_cast<size_t>(loaded->local_node_count));
      for (int local_index = 0; local_index < loaded->local_node_count;
           ++local_index) {
        const int global_node_id =
            loaded->local_to_global.empty()
                ? local_index
                : loaded->local_to_global[static_cast<size_t>(local_index)];
        result.full_labels.push_back(
            NodeLabel{global_node_id, local_index,
                      loaded->solver->getMinCutSolution(local_index)});
      }
    }
    markAlphaStateSolved(loaded);
    return result;
  }

  SolverStorageOptions storage_options_;
  std::unordered_map<int, LoadedPartition> partitions_;
};

class StreamingPartitionWorker final : public PartitionWorker {
public:
  struct Options {
    std::string storage_directory;
    std::uint64_t resident_byte_limit = 0;
    bool remove_storage_on_destroy = true;
    SolverStorageOptions solver_storage;
  };

  StreamingPartitionWorker() : StreamingPartitionWorker(Options{}) {}

  explicit StreamingPartitionWorker(Options options)
      : options_(std::move(options)) {
    if (options_.storage_directory.empty()) {
      auto name = std::string("mcpd3_streaming_worker_") +
                  std::to_string(
                      std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()) +
                  "_" +
                  std::to_string(reinterpret_cast<std::uintptr_t>(this));
      storage_directory_ = std::filesystem::temp_directory_path() / name;
      owns_storage_directory_ = true;
    } else {
      storage_directory_ = options_.storage_directory;
      owns_storage_directory_ = false;
    }
    std::filesystem::create_directories(storage_directory_);
  }

  ~StreamingPartitionWorker() override {
    evictAll();
    if (options_.remove_storage_on_destroy && owns_storage_directory_) {
      std::error_code ec;
      std::filesystem::remove_all(storage_directory_, ec);
    }
  }

  void loadPartition(const PartitionPackage &package) override {
    PartitionPackage copy = package;
    loadPartition(std::move(copy));
  }

  void loadPartition(PartitionPackage &&package) override {
    validatePartitionPackage(package);
    if (partitions_.find(package.partition_id) != partitions_.end()) {
      throw std::runtime_error("partition id " +
                               std::to_string(package.partition_id) +
                               " is already loaded");
    }

    StoredPartition stored;
    stored.partition_id = package.partition_id;
    stored.local_node_count = package.local_node_count;
    stored.local_arc_count = static_cast<int>(package.arcs.size() / 2);
    stored.path = storage_directory_ /
                  ("partition_" + std::to_string(package.partition_id) +
                   ".bin");
    stored.warm_state_path =
        storage_directory_ /
        ("partition_" + std::to_string(package.partition_id) + ".warm");
    stored.constraint_endpoints = std::move(package.constraint_endpoints);
    std::sort(stored.constraint_endpoints.begin(),
              stored.constraint_endpoints.end(),
              [](const ConstraintEndpointBinding &lhs,
                 const ConstraintEndpointBinding &rhs) {
                return lhs.constraint_id < rhs.constraint_id;
              });
    stored.resident_bytes = estimateResidentBytes(stored);

    writePackagePayload(stored.path, package);
    stored.local_to_global = std::move(package.local_to_global);
    partitions_.emplace(stored.partition_id, std::move(stored));
  }

  PartitionSolveResult solveRound(
      const PartitionSolveRequest &request) override {
    if (partitions_.empty()) {
      throw std::runtime_error("partition must be loaded before solveRound");
    }
    auto &stored = storedPartitionForRequest(request);
    const bool was_resident = static_cast<bool>(stored.resident_worker);
    applyAlphaUpdates(&stored, request.alpha_updates);
    auto *worker = materializePartition(&stored);

    PartitionSolveRequest forwarded = request;
    if (!was_resident) {
      forwarded.alpha_updates.clear();
    }
    auto result = worker->solveRound(forwarded);
    markAlphaStateSolved(&stored);
    stored.last_solution = worker->minCutSolution(stored.partition_id);
    stored.has_solution = true;
    stored.last_used = ++use_clock_;
    return result;
  }

  std::vector<PartitionSolveResult> solveRoundBatch(
      const std::vector<PartitionSolveRequest> &requests) override {
    std::vector<PartitionSolveResult> results;
    results.reserve(requests.size());
    std::unordered_set<int> seen_partition_ids;
    for (const auto &request : requests) {
      auto &stored = storedPartitionForRequest(request);
      if (!seen_partition_ids.insert(stored.partition_id).second) {
        throw std::runtime_error(
            "batch solve requests must target distinct partitions");
      }
    }
    for (const auto &request : requests) {
      results.push_back(solveRound(request));
    }
    return results;
  }

  void scaleObjective(long factor,
                      bool saturate_capacity_overflow = false) override {
    if (partitions_.empty()) {
      throw std::runtime_error("partition must be loaded before scaleObjective");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    for (auto &[partition_id, stored] : partitions_) {
      (void)partition_id;
      for (auto &binding : stored.constraint_endpoints) {
        binding.alpha = checkedScaleWorkerLagrange(binding.alpha, factor);
        binding.last_alpha =
            checkedScaleWorkerLagrange(binding.last_alpha, factor);
      }
      auto package = readPackagePayload(stored);
      for (auto &capacity : package.arc_capacities) {
        capacity = checkedScaleWorkerCapacity(
            capacity, factor, saturate_capacity_overflow);
      }
      for (auto &capacity : package.terminal_capacities) {
        capacity = checkedScaleWorkerCapacity(
            capacity, factor, saturate_capacity_overflow);
      }
      writePackagePayload(stored.path, package);
      if (stored.resident_worker) {
        stored.resident_worker->scaleObjective(factor,
                                               saturate_capacity_overflow);
      } else {
        invalidateWarmState(&stored);
      }
    }
  }

  void scaleObjectivePartitions(
      const std::vector<int> &partition_ids, long factor,
      bool saturate_capacity_overflow = false) override {
    if (partition_ids.empty()) {
      throw std::runtime_error(
          "objective scaling requires at least one partition id");
    }
    if (factor <= 0) {
      throw std::runtime_error("objective scale factor must be positive");
    }
    std::unordered_set<int> seen;
    std::vector<StoredPartition *> stored_partitions;
    stored_partitions.reserve(partition_ids.size());
    for (const int partition_id : partition_ids) {
      if (!seen.insert(partition_id).second) {
        throw std::runtime_error("duplicate objective scale partition id " +
                                 std::to_string(partition_id));
      }
      auto find_iter = partitions_.find(partition_id);
      if (find_iter == partitions_.end()) {
        throw std::runtime_error("unknown partition id " +
                                 std::to_string(partition_id));
      }
      stored_partitions.push_back(&find_iter->second);
    }
    for (auto *stored_ptr : stored_partitions) {
      auto &stored = *stored_ptr;
      const int partition_id = stored.partition_id;
      for (auto &binding : stored.constraint_endpoints) {
        binding.alpha = checkedScaleWorkerLagrange(binding.alpha, factor);
        binding.last_alpha =
            checkedScaleWorkerLagrange(binding.last_alpha, factor);
      }
      auto package = readPackagePayload(stored);
      for (auto &capacity : package.arc_capacities) {
        capacity = checkedScaleWorkerCapacity(
            capacity, factor, saturate_capacity_overflow);
      }
      for (auto &capacity : package.terminal_capacities) {
        capacity = checkedScaleWorkerCapacity(
            capacity, factor, saturate_capacity_overflow);
      }
      writePackagePayload(stored.path, package);
      if (stored.resident_worker) {
        stored.resident_worker->scaleObjectivePartitions(
            {partition_id}, factor, saturate_capacity_overflow);
      } else {
        invalidateWarmState(&stored);
      }
    }
  }

  void replacePartitionCapacities(
      const PartitionCapacityUpdate &update) override {
    replacePartitionCapacitiesFor(update.partition_id, update);
  }

  void replacePartitionCapacitiesFor(
      int target_partition_id,
      const PartitionCapacityUpdate &update) override {
    auto find_iter = partitions_.find(target_partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown partition id " +
                               std::to_string(target_partition_id));
    }
    auto &stored = find_iter->second;
    auto *worker = materializePartition(&stored);
    worker->replacePartitionCapacitiesFor(target_partition_id, update);
    auto package = readPackagePayload(stored);
    package.arc_capacities = update.arc_capacities;
    package.terminal_capacities = update.terminal_capacities;
    writePackagePayload(stored.path, package);
    stored.last_solution = worker->minCutSolution(stored.partition_id);
    stored.has_solution = true;
    stored.has_warm_state = false;
    std::error_code error;
    std::filesystem::remove(stored.warm_state_path, error);
  }

  std::size_t fullLabelCount(int partition_id) const override {
    auto find_iter = partitions_.find(partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown partition id " +
                               std::to_string(partition_id));
    }
    return static_cast<std::size_t>(find_iter->second.local_node_count);
  }

  void copyFullLabels(int partition_id, std::size_t offset,
                      NodeLabel *destination, std::size_t count) override {
    auto find_iter = partitions_.find(partition_id);
    if (find_iter == partitions_.end()) {
      throw std::runtime_error("unknown partition id " +
                               std::to_string(partition_id));
    }
    auto &stored = find_iter->second;
    auto *worker = materializePartition(&stored);
    worker->copyFullLabels(partition_id, offset, destination, count);
    stored.last_used = ++use_clock_;
  }

  std::uint64_t residentBytesForTesting() const { return resident_bytes_; }
  long warmStateWriteCountForTesting() const { return warm_state_write_count_; }
  long warmStateRestoreCountForTesting() const {
    return warm_state_restore_count_;
  }
  long residentPartitionCountForTesting() const {
    long count = 0;
    for (const auto &entry : partitions_) {
      if (entry.second.resident_worker) {
        ++count;
      }
    }
    return count;
  }

private:
  struct StoredPartition {
    int partition_id = -1;
    int local_node_count = 0;
    int local_arc_count = 0;
    std::filesystem::path path;
    std::filesystem::path warm_state_path;
    std::uint64_t resident_bytes = 0;
    std::uint64_t last_used = 0;
    bool has_solution = false;
    bool has_warm_state = false;
    SolverArray<int> local_to_global;
    std::vector<ConstraintEndpointBinding> constraint_endpoints;
    std::vector<int> last_solution;
    std::unique_ptr<InProcessPartitionWorker> resident_worker;
  };

  static std::uint64_t endpointBytes(const StoredPartition &stored) {
    return static_cast<std::uint64_t>(stored.constraint_endpoints.size()) *
           static_cast<std::uint64_t>(sizeof(ConstraintEndpointBinding));
  }

  static std::uint64_t localToGlobalBytes(const StoredPartition &stored) {
    return static_cast<std::uint64_t>(stored.local_to_global.size()) *
           static_cast<std::uint64_t>(sizeof(int));
  }

  static std::uint64_t estimateResidentBytes(const StoredPartition &stored) {
    const auto estimate = PrimalDualMinCutSolver::estimateMemoryBytes(
        stored.local_node_count, stored.local_arc_count);
    return static_cast<std::uint64_t>(estimate.total_bytes) +
           endpointBytes(stored) + localToGlobalBytes(stored);
  }

  template <typename T>
  static void writeScalar(std::ostream &out, const T &value,
                          const std::string &name) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(T));
    if (!out) {
      throw std::runtime_error("failed to write " + name);
    }
  }

  template <typename T>
  static T readScalar(std::istream &in, const std::string &name) {
    T value{};
    in.read(reinterpret_cast<char *>(&value), sizeof(T));
    if (!in) {
      throw std::runtime_error("failed to read " + name);
    }
    return value;
  }

  template <typename Container>
  static void writeVector(std::ostream &out, const Container &values,
                          const std::string &name) {
    using T = typename Container::value_type;
    const std::uint64_t size = values.size();
    writeScalar(out, size, name + " size");
    if (!values.empty()) {
      out.write(reinterpret_cast<const char *>(values.data()),
                static_cast<std::streamsize>(values.size() * sizeof(T)));
      if (!out) {
        throw std::runtime_error("failed to write " + name);
      }
    }
  }

  template <typename T>
  static std::vector<T> readVector(std::istream &in,
                                   const std::string &name) {
    const auto size = readScalar<std::uint64_t>(in, name + " size");
    if (size >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error(name + " is too large");
    }
    std::vector<T> values(static_cast<std::size_t>(size));
    if (!values.empty()) {
      in.read(reinterpret_cast<char *>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(T)));
      if (!in) {
        throw std::runtime_error("failed to read " + name);
      }
    }
    return values;
  }

  template <typename Integer>
  static void writeInteger(std::ostream &out, const Integer &value,
                           const std::string &name) {
    if constexpr (std::is_trivially_copyable_v<Integer>) {
      writeScalar(out, value, name);
    } else {
      const std::string text = integer_to_string(value);
      const std::uint64_t size = text.size();
      writeScalar(out, size, name + " size");
      out.write(text.data(), static_cast<std::streamsize>(text.size()));
      if (!out) {
        throw std::runtime_error("failed to write " + name);
      }
    }
  }

  template <typename Integer, typename Parser>
  static Integer readInteger(std::istream &in, const std::string &name,
                             Parser parser) {
    if constexpr (std::is_trivially_copyable_v<Integer>) {
      return readScalar<Integer>(in, name);
    } else {
      const auto size = readScalar<std::uint64_t>(in, name + " size");
      if (size > static_cast<std::uint64_t>(
                     std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error(name + " is too large");
      }
      std::string text(static_cast<std::size_t>(size), '\0');
      in.read(text.data(), static_cast<std::streamsize>(text.size()));
      if (!in) {
        throw std::runtime_error("failed to read " + name);
      }
      return parser(text);
    }
  }

  template <typename Container>
  static void writeIntegerVector(std::ostream &out,
                                 const Container &values,
                                 const std::string &name) {
    using Integer = typename Container::value_type;
    const std::uint64_t size = values.size();
    writeScalar(out, size, name + " size");
    if constexpr (std::is_trivially_copyable_v<Integer>) {
      if (!values.empty()) {
        out.write(reinterpret_cast<const char *>(values.data()),
                  static_cast<std::streamsize>(values.size() *
                                               sizeof(Integer)));
        if (!out) {
          throw std::runtime_error("failed to write " + name);
        }
      }
    } else {
      for (const auto &value : values) {
        writeInteger(out, value, name + " value");
      }
    }
  }

  template <typename Integer, typename Parser>
  static std::vector<Integer>
  readIntegerVector(std::istream &in, const std::string &name, Parser parser) {
    const auto size = readScalar<std::uint64_t>(in, name + " size");
    if (size >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error(name + " is too large");
    }
    std::vector<Integer> values;
    values.reserve(static_cast<std::size_t>(size));
    if constexpr (std::is_trivially_copyable_v<Integer>) {
      values.resize(static_cast<std::size_t>(size));
      if (!values.empty()) {
        in.read(reinterpret_cast<char *>(values.data()),
                static_cast<std::streamsize>(values.size() *
                                             sizeof(Integer)));
        if (!in) {
          throw std::runtime_error("failed to read " + name);
        }
      }
    } else {
      for (std::uint64_t index = 0; index < size; ++index) {
        values.push_back(readInteger<Integer>(in, name + " value", parser));
      }
    }
    return values;
  }

  template <typename Container>
  static void writeIntVector(std::ostream &out, const Container &values,
                             const std::string &name) {
    writeVector(out, values, name);
  }

  static std::vector<int> readIntVector(std::istream &in,
                                        const std::string &name) {
    return readVector<int>(in, name);
  }

  static void writePackagePayload(const std::filesystem::path &path,
                                  const PartitionPackage &package) {
    const auto tmp_path = path.string() + ".tmp";
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("failed to open streaming package file for " +
                               path.string());
    }
    const std::uint32_t magic = 0x4d435033;
    const std::uint32_t version = 4;
    writeScalar(out, magic, "package magic");
    writeScalar(out, version, "package version");
    writeScalar(out, package.partition_id, "partition id");
    writeScalar(out, package.local_node_count, "local node count");
    writeIntVector(out, package.arcs, "arcs");
    writeIntegerVector(out, package.arc_capacities, "arc capacities");
    writeIntegerVector(out, package.terminal_capacities,
                       "terminal capacities");
    writeIntVector(out, package.local_to_global, "local to global");
    writeScalar(out, package.objective_multiplier, "objective multiplier");
    writeScalar(out,
                static_cast<std::uint32_t>(package.canonical_cut_selection),
                "canonical cut selection");
    writeScalar(out,
                static_cast<std::uint8_t>(
                    package.force_full_mincut_recompute ? 1 : 0),
                "force full mincut recompute");
    writeIntVector(out, package.reference_cut_labels,
                   "reference cut labels");
    writeScalar(out,
                static_cast<std::uint32_t>(package.reference_cut_selection),
                "reference cut selection");
    writeScalar(out, package.reference_cut_check_interval,
                "reference cut check interval");
    out.close();
    if (!out) {
      throw std::runtime_error("failed to flush streaming package file " +
                               path.string());
    }
    std::filesystem::rename(tmp_path, path);
  }

  static void writeWarmState(
      const std::filesystem::path &path,
      const PrimalDualMinCutSolver::WarmState &state) {
    const auto tmp_path = path.string() + ".tmp";
    std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
    if (!out) {
      throw std::runtime_error("failed to open streaming warm-state file for " +
                               path.string());
    }
    const std::uint32_t magic = 0x4d435357;
    const std::uint32_t version = 4;
    writeScalar(out, magic, "warm state magic");
    writeScalar(out, version, "warm state version");
    writeIntegerVector(out, state.v_flow, "v flow");
    writeIntegerVector(out, state.d_flow, "d flow");
    writeVector(out, state.x, "min cut labels");
    writeScalar(out, static_cast<std::uint8_t>(state.is_first_iteration ? 1 : 0),
                "is first iteration");
    writeScalar(out,
                static_cast<std::uint8_t>(
                    state.is_first_iteration_of_new_scale ? 1 : 0),
                "is first iteration of new scale");
    writeScalar(out, static_cast<std::uint8_t>(state.has_solution ? 1 : 0),
                "has solution");
    writeInteger(out, state.mincut_value, "mincut value");
    writeIntegerVector(out, state.cached_lagrange_multipliers,
                       "cached lagrange multipliers");
    writeIntegerVector(out, state.cached_last_lagrange_multipliers,
                       "cached last lagrange multipliers");
    writeInteger(out, state.regularization_str, "regularization strength");
    writeInteger(out, state.last_regularization_budget,
                 "last regularization budget");
    writeInteger(out, state.last_regularization_contribution,
                 "last regularization contribution");
    writeScalar(out, state.last_regularization_anchor_sink_count,
                "last regularization anchor sink count");
    writeScalar(out, state.last_regularization_active_sink_count,
                "last regularization active sink count");
    writeIntegerVector(out, state.regularization_weights,
                       "regularization weights");
    const auto &graph_state = state.maxflow_graph_state;
    writeScalar(out, graph_state.node_num, "warm graph node count");
    writeScalar(out, graph_state.arc_num, "warm graph arc count");
    writeInteger(out, graph_state.flow, "warm graph flow");
    writeScalar(out, graph_state.maxflow_iteration,
                "warm graph maxflow iteration");
    writeScalar(out, graph_state.time, "warm graph time");
    writeIntegerVector(out, graph_state.node_tr_caps,
                       "warm graph node tr caps");
    writeVector(out, graph_state.node_parent_arc_indices,
                "warm graph node parent arc indices");
    writeVector(out, graph_state.node_timestamps,
                "warm graph node timestamps");
    writeVector(out, graph_state.node_distances, "warm graph node distances");
    writeVector(out, graph_state.node_is_sink, "warm graph node is sink");
    writeIntegerVector(out, graph_state.arc_residual_capacities,
                       "warm graph arc residual capacities");
    out.close();
    if (!out) {
      throw std::runtime_error("failed to flush streaming warm-state file " +
                               path.string());
    }
    std::filesystem::rename(tmp_path, path);
  }

  static PrimalDualMinCutSolver::WarmState
  readWarmState(const std::filesystem::path &path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
      throw std::runtime_error("failed to open streaming warm-state file " +
                               path.string());
    }
    const auto magic = readScalar<std::uint32_t>(in, "warm state magic");
    const auto version = readScalar<std::uint32_t>(in, "warm state version");
    if (magic != 0x4d435357 || (version != 3 && version != 4)) {
      throw std::runtime_error("invalid streaming warm-state file " +
                               path.string());
    }
    PrimalDualMinCutSolver::WarmState state;
    state.v_flow = readIntegerVector<Capacity>(in, "v flow", parse_capacity);
    state.d_flow =
        readIntegerVector<NodeFlow>(in, "d flow", parse_node_flow);
    state.x = readVector<int>(in, "min cut labels");
    state.is_first_iteration =
        readScalar<std::uint8_t>(in, "is first iteration") != 0;
    state.is_first_iteration_of_new_scale =
        readScalar<std::uint8_t>(in, "is first iteration of new scale") != 0;
    state.has_solution = readScalar<std::uint8_t>(in, "has solution") != 0;
    state.mincut_value =
        readInteger<Objective>(in, "mincut value", parse_objective);
    state.cached_lagrange_multipliers =
        readIntegerVector<Lagrange>(in, "cached lagrange multipliers",
                                    parse_lagrange);
    state.cached_last_lagrange_multipliers =
        readIntegerVector<Lagrange>(in, "cached last lagrange multipliers",
                                    parse_lagrange);
    state.regularization_str =
        readInteger<Capacity>(in, "regularization strength", parse_capacity);
    state.last_regularization_budget =
        readInteger<Objective>(in, "last regularization budget",
                               parse_objective);
    state.last_regularization_contribution =
        readInteger<Objective>(in, "last regularization contribution",
                               parse_objective);
    state.last_regularization_anchor_sink_count =
        readScalar<long>(in, "last regularization anchor sink count");
    state.last_regularization_active_sink_count =
        readScalar<long>(in, "last regularization active sink count");
    if (version == 3) {
      const auto anchors =
          readVector<unsigned char>(in, "regularization anchor sink");
      state.regularization_weights.reserve(anchors.size());
      for (const unsigned char anchor : anchors) {
        state.regularization_weights.push_back(
            anchor ? widen_capacity(state.regularization_str) : Objective{0});
      }
    } else {
      state.regularization_weights = readIntegerVector<Objective>(
          in, "regularization weights", parse_objective);
    }
    auto &graph_state = state.maxflow_graph_state;
    graph_state.node_num = readScalar<int>(in, "warm graph node count");
    graph_state.arc_num = readScalar<int>(in, "warm graph arc count");
    graph_state.flow =
        readInteger<Objective>(in, "warm graph flow", parse_objective);
    graph_state.maxflow_iteration =
        readScalar<int>(in, "warm graph maxflow iteration");
    graph_state.time = readScalar<long>(in, "warm graph time");
    graph_state.node_tr_caps = readIntegerVector<TerminalResidual>(
        in, "warm graph node tr caps", parse_terminal_residual);
    graph_state.node_parent_arc_indices =
        readVector<int>(in, "warm graph node parent arc indices");
    graph_state.node_timestamps =
        readVector<long>(in, "warm graph node timestamps");
    graph_state.node_distances =
        readVector<int>(in, "warm graph node distances");
    graph_state.node_is_sink =
        readVector<unsigned char>(in, "warm graph node is sink");
    graph_state.arc_residual_capacities = readIntegerVector<Capacity>(
        in, "warm graph arc residual capacities", parse_capacity);
    return state;
  }

  static PartitionPackage readPackagePayload(const StoredPartition &stored) {
    std::ifstream in(stored.path, std::ios::binary);
    if (!in) {
      throw std::runtime_error("failed to open streaming package file " +
                               stored.path.string());
    }
    const auto magic = readScalar<std::uint32_t>(in, "package magic");
    const auto version = readScalar<std::uint32_t>(in, "package version");
    if (magic != 0x4d435033 || (version != 3 && version != 4)) {
      throw std::runtime_error("invalid streaming package file " +
                               stored.path.string());
    }

    PartitionPackage package;
    package.partition_id = readScalar<int>(in, "partition id");
    package.local_node_count = readScalar<int>(in, "local node count");
    package.arcs = readIntVector(in, "arcs");
    package.arc_capacities = readIntegerVector<Capacity>(
        in, "arc capacities", parse_capacity);
    package.terminal_capacities = readIntegerVector<Capacity>(
        in, "terminal capacities", parse_capacity);
    package.local_to_global = readIntVector(in, "local to global");
    if (version >= 4) {
      package.objective_multiplier =
          readScalar<long>(in, "objective multiplier");
      const auto canonical_selection =
          readScalar<std::uint32_t>(in, "canonical cut selection");
      if (canonical_selection > static_cast<std::uint32_t>(
                                    CanonicalCutSelection::MAXIMUM_LABELS)) {
        throw std::runtime_error(
            "invalid streaming canonical cut selection");
      }
      package.canonical_cut_selection =
          static_cast<CanonicalCutSelection>(canonical_selection);
      package.force_full_mincut_recompute =
          readScalar<std::uint8_t>(in, "force full mincut recompute") != 0;
      package.reference_cut_labels =
          readIntVector(in, "reference cut labels");
      const auto reference_selection =
          readScalar<std::uint32_t>(in, "reference cut selection");
      if (reference_selection >
          static_cast<std::uint32_t>(
              ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL)) {
        throw std::runtime_error(
            "invalid streaming reference cut selection");
      }
      package.reference_cut_selection =
          static_cast<ReferenceCutSelection>(reference_selection);
      package.reference_cut_check_interval =
          readScalar<long>(in, "reference cut check interval");
    }
    package.constraint_endpoints = stored.constraint_endpoints;
    validatePartitionPackage(package);
    return package;
  }

  StoredPartition &storedPartitionForRequest(
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

  ConstraintEndpointBinding &endpointForUpdate(StoredPartition *stored,
                                               const AlphaUpdate &update) {
    auto iter = std::lower_bound(
        stored->constraint_endpoints.begin(), stored->constraint_endpoints.end(),
        update.constraint_id,
        [](const ConstraintEndpointBinding &binding, int constraint_id) {
          return binding.constraint_id < constraint_id;
        });
    if (iter == stored->constraint_endpoints.end() ||
        iter->constraint_id != update.constraint_id) {
      throw std::runtime_error("unknown alpha update constraint id " +
                               std::to_string(update.constraint_id));
    }
    return *iter;
  }

  void applyAlphaUpdates(StoredPartition *stored,
                         const std::vector<AlphaUpdate> &updates) {
    for (const auto &update : updates) {
      auto &binding = endpointForUpdate(stored, update);
      binding.last_alpha = binding.alpha;
      binding.alpha = update.alpha;
    }
  }

  void markAlphaStateSolved(StoredPartition *stored) {
    for (auto &binding : stored->constraint_endpoints) {
      binding.last_alpha = binding.alpha;
    }
  }

  InProcessPartitionWorker *materializePartition(StoredPartition *stored) {
    if (stored->resident_worker) {
      return stored->resident_worker.get();
    }
    evictUntilFits(stored->resident_bytes);
    auto package = readPackagePayload(*stored);
    stored->resident_worker =
        std::make_unique<InProcessPartitionWorker>(options_.solver_storage);
    stored->resident_worker->loadPartition(std::move(package));
    if (stored->has_warm_state) {
      stored->resident_worker->restoreWarmState(
          stored->partition_id, readWarmState(stored->warm_state_path));
      ++warm_state_restore_count_;
    } else if (stored->has_solution) {
      stored->resident_worker->restoreMinCutSolution(stored->partition_id,
                                                     stored->last_solution);
    }
    resident_bytes_ += stored->resident_bytes;
    stored->last_used = ++use_clock_;
    return stored->resident_worker.get();
  }

  void evictUntilFits(std::uint64_t incoming_bytes) {
    if (options_.resident_byte_limit == 0) {
      return;
    }
    while (resident_bytes_ + incoming_bytes > options_.resident_byte_limit) {
      auto evict_iter = partitions_.end();
      for (auto iter = partitions_.begin(); iter != partitions_.end(); ++iter) {
        if (!iter->second.resident_worker) {
          continue;
        }
        if (evict_iter == partitions_.end() ||
            iter->second.last_used < evict_iter->second.last_used) {
          evict_iter = iter;
        }
      }
      if (evict_iter == partitions_.end()) {
        return;
      }
      evictResident(&evict_iter->second);
    }
  }

  void evictResident(StoredPartition *stored, bool persist_warm_state = true) {
    if (!stored->resident_worker) {
      return;
    }
    if (persist_warm_state) {
      writeWarmState(stored->warm_state_path,
                     stored->resident_worker->warmState(stored->partition_id));
      stored->has_warm_state = true;
      ++warm_state_write_count_;
    }
    stored->resident_worker.reset();
    resident_bytes_ -= stored->resident_bytes;
  }

  void invalidateWarmState(StoredPartition *stored) {
    stored->has_warm_state = false;
    std::error_code ec;
    std::filesystem::remove(stored->warm_state_path, ec);
  }

  void evictAll() {
    for (auto &[partition_id, stored] : partitions_) {
      (void)partition_id;
      evictResident(&stored, /*persist_warm_state=*/false);
    }
  }

  Options options_;
  std::filesystem::path storage_directory_;
  bool owns_storage_directory_ = false;
  std::uint64_t resident_bytes_ = 0;
  std::uint64_t use_clock_ = 0;
  long warm_state_write_count_ = 0;
  long warm_state_restore_count_ = 0;
  std::unordered_map<int, StoredPartition> partitions_;
};

} // namespace mcpd3
