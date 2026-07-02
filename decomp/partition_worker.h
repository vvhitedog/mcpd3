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

struct PartitionWorkerResourceEstimate {
  int cpu_count = 1;
  long ram_gb = 0;
};

inline long checkedScaleWorkerLong(long value, long scale) {
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

inline int checkedScaleWorkerInt(int value, long scale,
                                 bool saturate_capacity_overflow = false) {
  const long result = checkedScaleWorkerLong(value, scale);
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

inline void validatePartitionPackage(const PartitionPackage &package) {
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
};

class InProcessPartitionWorker final : public PartitionWorker {
public:
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
    loaded.constraint_endpoints = std::move(package.constraint_endpoints);
    loaded.solver = std::make_unique<PrimalDualMinCutSolver>(
        local_node_count, arc_count, std::move(package.arcs),
        std::move(package.arc_capacities),
        std::move(package.terminal_capacities));

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
        constraint_arc.alpha = checkedScaleLong(constraint_arc.alpha, factor);
        constraint_arc.last_alpha =
            checkedScaleLong(constraint_arc.last_alpha, factor);
      }
      loaded.solver->scaleProblem(factor, saturate_capacity_overflow);
    }
  }

private:
  struct LoadedPartition {
    int partition_id = -1;
    int local_node_count = 0;
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
    markAlphaStateSolved(loaded);
    return result;
  }

  std::unordered_map<int, LoadedPartition> partitions_;
};

class StreamingPartitionWorker final : public PartitionWorker {
public:
  struct Options {
    std::string storage_directory;
    std::uint64_t resident_byte_limit = 0;
    bool remove_storage_on_destroy = true;
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
    stored.constraint_endpoints = std::move(package.constraint_endpoints);
    std::sort(stored.constraint_endpoints.begin(),
              stored.constraint_endpoints.end(),
              [](const ConstraintEndpointBinding &lhs,
                 const ConstraintEndpointBinding &rhs) {
                return lhs.constraint_id < rhs.constraint_id;
              });
    stored.resident_bytes = estimateResidentBytes(stored);

    writePackagePayload(stored.path, package);
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
        binding.alpha = checkedScaleWorkerLong(binding.alpha, factor);
        binding.last_alpha = checkedScaleWorkerLong(binding.last_alpha, factor);
      }
      auto package = readPackagePayload(stored);
      for (auto &capacity : package.arc_capacities) {
        capacity =
            checkedScaleWorkerInt(capacity, factor, saturate_capacity_overflow);
      }
      for (auto &capacity : package.terminal_capacities) {
        capacity =
            checkedScaleWorkerInt(capacity, factor, saturate_capacity_overflow);
      }
      writePackagePayload(stored.path, package);
      if (stored.resident_worker) {
        stored.resident_worker->scaleObjective(factor,
                                               saturate_capacity_overflow);
      }
    }
  }

  std::uint64_t residentBytesForTesting() const { return resident_bytes_; }
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
    std::uint64_t resident_bytes = 0;
    std::uint64_t last_used = 0;
    bool has_solution = false;
    std::vector<ConstraintEndpointBinding> constraint_endpoints;
    std::vector<int> last_solution;
    std::unique_ptr<InProcessPartitionWorker> resident_worker;
  };

  static std::uint64_t endpointBytes(const StoredPartition &stored) {
    return static_cast<std::uint64_t>(stored.constraint_endpoints.size()) *
           static_cast<std::uint64_t>(sizeof(ConstraintEndpointBinding));
  }

  static std::uint64_t estimateResidentBytes(const StoredPartition &stored) {
    const auto estimate = PrimalDualMinCutSolver::estimateMemoryBytes(
        stored.local_node_count, stored.local_arc_count);
    return static_cast<std::uint64_t>(estimate.total_bytes) +
           endpointBytes(stored);
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

  static void writeIntVector(std::ostream &out,
                             const std::vector<int> &values,
                             const std::string &name) {
    const std::uint64_t size = values.size();
    writeScalar(out, size, name + " size");
    if (!values.empty()) {
      out.write(reinterpret_cast<const char *>(values.data()),
                static_cast<std::streamsize>(values.size() * sizeof(int)));
      if (!out) {
        throw std::runtime_error("failed to write " + name);
      }
    }
  }

  static std::vector<int> readIntVector(std::istream &in,
                                        const std::string &name) {
    const auto size = readScalar<std::uint64_t>(in, name + " size");
    if (size >
        static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error(name + " is too large");
    }
    std::vector<int> values(static_cast<std::size_t>(size));
    if (!values.empty()) {
      in.read(reinterpret_cast<char *>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(int)));
      if (!in) {
        throw std::runtime_error("failed to read " + name);
      }
    }
    return values;
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
    const std::uint32_t version = 1;
    writeScalar(out, magic, "package magic");
    writeScalar(out, version, "package version");
    writeScalar(out, package.partition_id, "partition id");
    writeScalar(out, package.local_node_count, "local node count");
    writeIntVector(out, package.arcs, "arcs");
    writeIntVector(out, package.arc_capacities, "arc capacities");
    writeIntVector(out, package.terminal_capacities, "terminal capacities");
    out.close();
    if (!out) {
      throw std::runtime_error("failed to flush streaming package file " +
                               path.string());
    }
    std::filesystem::rename(tmp_path, path);
  }

  static PartitionPackage readPackagePayload(const StoredPartition &stored) {
    std::ifstream in(stored.path, std::ios::binary);
    if (!in) {
      throw std::runtime_error("failed to open streaming package file " +
                               stored.path.string());
    }
    const auto magic = readScalar<std::uint32_t>(in, "package magic");
    const auto version = readScalar<std::uint32_t>(in, "package version");
    if (magic != 0x4d435033 || version != 1) {
      throw std::runtime_error("invalid streaming package file " +
                               stored.path.string());
    }

    PartitionPackage package;
    package.partition_id = readScalar<int>(in, "partition id");
    package.local_node_count = readScalar<int>(in, "local node count");
    package.arcs = readIntVector(in, "arcs");
    package.arc_capacities = readIntVector(in, "arc capacities");
    package.terminal_capacities = readIntVector(in, "terminal capacities");
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
    stored->resident_worker = std::make_unique<InProcessPartitionWorker>();
    stored->resident_worker->loadPartition(std::move(package));
    if (stored->has_solution) {
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

  void evictResident(StoredPartition *stored) {
    if (!stored->resident_worker) {
      return;
    }
    stored->resident_worker.reset();
    resident_bytes_ -= stored->resident_bytes;
  }

  void evictAll() {
    for (auto &[partition_id, stored] : partitions_) {
      (void)partition_id;
      evictResident(&stored);
    }
  }

  Options options_;
  std::filesystem::path storage_directory_;
  bool owns_storage_directory_ = false;
  std::uint64_t resident_bytes_ = 0;
  std::uint64_t use_clock_ = 0;
  std::unordered_map<int, StoredPartition> partitions_;
};

} // namespace mcpd3
