// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#include <cstdlib>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <limits>
#include <new>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <unistd.h>

#include <decomp/dualdecomp.h>
#include <decomp/halo_partition.h>
#include <decomp/partition_coordinator.h>
#include <decomp/partition_worker.h>
#include <decomp/speculative_cycle_replay.h>
#include <graph/dimacs.h>
#include <primaldual/mcpd3.h>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

std::string captureStderr(const std::function<void()> &function) {
  std::fflush(stderr);
  const int saved_stderr = ::dup(STDERR_FILENO);
  require(saved_stderr >= 0, "failed to duplicate stderr");
  FILE *capture = std::tmpfile();
  require(capture != nullptr, "failed to create stderr capture file");
  require(::dup2(::fileno(capture), STDERR_FILENO) >= 0,
          "failed to redirect stderr");

  try {
    function();
  } catch (...) {
    std::fflush(stderr);
    (void)::dup2(saved_stderr, STDERR_FILENO);
    ::close(saved_stderr);
    std::fclose(capture);
    throw;
  }

  std::fflush(stderr);
  require(std::fseek(capture, 0, SEEK_END) == 0,
          "failed to seek captured stderr");
  const long capture_size = std::ftell(capture);
  require(capture_size >= 0, "failed to size captured stderr");
  std::rewind(capture);
  std::string output(static_cast<size_t>(capture_size), '\0');
  require(output.empty() ||
              std::fread(output.data(), 1, output.size(), capture) ==
                  output.size(),
          "failed to read captured stderr");
  require(::dup2(saved_stderr, STDERR_FILENO) >= 0,
          "failed to restore stderr");
  ::close(saved_stderr);
  std::fclose(capture);
  return output;
}

template <typename Function>
void requireThrows(Function function, const std::string &message) {
  bool threw = false;
  try {
    function();
  } catch (const std::exception &) {
    threw = true;
  }
  require(threw, message);
}

void haloDepthOneRetainsLegacyOwnerLayout() {
  const auto layout = mcpd3::buildHaloPartitionLayout(
      /*partition_count=*/2, /*node_count=*/6,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3, 3, 4, 4, 5},
      /*partition_labels=*/std::vector<int>{0, 0, 0, 1, 1, 1},
      /*halo_depth=*/1);

  require(layout.objective_multiplier == 1,
          "h1 compatibility layout must not rescale the objective");
  require(layout.arc_partitions ==
              std::vector<std::vector<int>>({{0}, {0}, {0}, {1}, {1}}),
          "h1 compatibility layout must retain one-owner edge placement");
  require(layout.node_partitions ==
              std::vector<std::vector<int>>(
                  {{0}, {0}, {0}, {0, 1}, {1}, {1}}),
          "h1 compatibility layout must retain legacy boundary clones");
}

void finiteHaloLayoutUsesBfsMembershipAndExactMultiplicity() {
  const auto layout = mcpd3::buildHaloPartitionLayout(
      /*partition_count=*/2, /*node_count=*/6,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3, 3, 4, 4, 5},
      /*partition_labels=*/std::vector<int>{0, 0, 0, 1, 1, 1},
      /*halo_depth=*/2);

  require(layout.objective_multiplier == 2,
          "two overlapping h2 halos should require objective multiplier 2");
  require(layout.node_partitions ==
              std::vector<std::vector<int>>(
                  {{0}, {0, 1}, {0, 1}, {0, 1}, {0, 1}, {1}}),
          "h2 node memberships must match two-step BFS from each core");
  require(layout.arc_partitions ==
              std::vector<std::vector<int>>(
                  {{0}, {0, 1}, {0, 1}, {0, 1}, {1}}),
          "an h2 edge must appear in every halo containing both endpoints");
}

void haloObjectiveMultiplierUsesNodeAndEdgeLcm() {
  const auto layout = mcpd3::buildHaloPartitionLayout(
      /*partition_count=*/3, /*node_count=*/6,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3, 3, 4, 4, 5},
      /*partition_labels=*/std::vector<int>{0, 0, 1, 1, 2, 2},
      /*halo_depth=*/2);

  require(layout.objective_multiplier == 6,
          "mixed two-way and three-way halo terms require LCM 6");
  require(layout.node_partitions[2] == std::vector<int>({0, 1, 2}) &&
              layout.arc_partitions[2] == std::vector<int>({0, 1, 2}),
          "center node and edge should occur in all three h2 halos");
  require(layout.node_partitions[0] == std::vector<int>({0, 1}) &&
              layout.arc_partitions[0] == std::vector<int>({0, 1}),
          "outer node and edge should occur in exactly two h2 halos");
}

void infiniteHaloContainsEveryNodeAndEdge() {
  const auto layout = mcpd3::buildHaloPartitionLayout(
      /*partition_count=*/2, /*node_count=*/4,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3},
      /*partition_labels=*/std::vector<int>{0, 0, 1, 1},
      mcpd3::kInfiniteHaloDepth);

  require(layout.objective_multiplier == 2,
          "two full-graph halos should require multiplier 2");
  for (const auto &memberships : layout.node_partitions) {
    require(memberships == std::vector<int>({0, 1}),
            "infinite halo must contain every node in every partition");
  }
  for (const auto &memberships : layout.arc_partitions) {
    require(memberships == std::vector<int>({0, 1}),
            "infinite halo must contain every edge in every partition");
  }
}

void haloLayoutRejectsInvalidInputs() {
  requireThrows(
      [] {
        (void)mcpd3::buildHaloPartitionLayout(
            2, 2, std::vector<int>{0, 1}, std::vector<int>{0, 1}, 0);
      },
      "zero halo depth should be rejected");
  requireThrows(
      [] {
        (void)mcpd3::buildHaloPartitionLayout(
            2, 2, std::vector<int>{0, 2}, std::vector<int>{0, 1}, 2);
      },
      "out-of-range halo arc endpoints should be rejected");
  requireThrows(
      [] {
        (void)mcpd3::buildHaloPartitionLayout(
            2, 2, std::vector<int>{0, 1}, std::vector<int>{0}, 2);
      },
      "wrong-sized halo partition labels should be rejected");
  requireThrows(
      [] {
        (void)mcpd3::checkedHaloLcm(std::numeric_limits<long>::max(), 2);
      },
      "halo objective multiplier overflow should be rejected");
}

void lowerBoundCertificateSubtractsOnlyRegularizationSlack() {
  require(mcpd3::regularizedObjectiveRaw(/*original_objective_raw=*/100,
                                         /*regularization_contribution_raw=*/7) ==
              107,
          "regularized objective should add the actual contribution");
  require(mcpd3::certifiedOriginalLowerBoundRaw(
              /*original_objective_raw=*/100,
              /*regularization_contribution_raw=*/7,
              /*regularization_budget_raw=*/10) == 97,
          "certified lower bound should subtract only budget slack");
  require(mcpd3::certifiedOriginalLowerBoundRaw(
              /*original_objective_raw=*/-20,
              /*regularization_contribution_raw=*/3,
              /*regularization_budget_raw=*/5) == -22,
          "certificate arithmetic should support negative local objectives");

  bool add_threw = false;
  try {
    (void)mcpd3::regularizedObjectiveRaw(
        std::numeric_limits<mcpd3::Objective>::max(), mcpd3::Objective(1));
  } catch (const std::overflow_error &) {
    add_threw = true;
  }
  require(add_threw ==
              (mcpd3::integer_is_bounded<mcpd3::Objective>() &&
               !mcpd3::legacy_32bit_dd_replay_enabled()),
          "regularized objective overflow behavior should match precision");

  bool subtract_threw = false;
  try {
    (void)mcpd3::certifiedOriginalLowerBoundRaw(
        std::numeric_limits<mcpd3::Objective>::min(), mcpd3::Objective(0),
        mcpd3::Objective(1));
  } catch (const std::overflow_error &) {
    subtract_threw = true;
  }
  require(subtract_threw ==
              (mcpd3::integer_is_bounded<mcpd3::Objective>() &&
               !mcpd3::legacy_32bit_dd_replay_enabled()),
          "certificate underflow behavior should match precision");
}

void solverMemoryEstimateReportsBkAndVectorBytes() {
  using GraphType =
      Graph<mcpd3::Capacity, mcpd3::TerminalResidual, mcpd3::Objective>;
  require(GraphType::estimated_node_array_bytes(1) ==
              GraphType::estimated_node_array_bytes(16),
          "BK node estimate should include constructor minimum capacity");
  require(GraphType::estimated_arc_array_bytes(1) ==
              GraphType::estimated_arc_array_bytes(16),
          "BK arc estimate should include constructor minimum capacity");

  const auto estimate =
      mcpd3::PrimalDualMinCutSolver::estimateMemoryBytes(/*nnode=*/2,
                                                         /*narc=*/1);
  require(estimate.bk_node_bytes > 0, "BK node estimate should be positive");
  require(estimate.bk_arc_bytes > 0, "BK arc estimate should be positive");
  require(estimate.bk_total_bytes ==
              estimate.bk_node_bytes + estimate.bk_arc_bytes,
          "BK total estimate should sum node and arc arrays");
  require(estimate.solver_vector_bytes ==
              7 * sizeof(int) + 5 * sizeof(mcpd3::Capacity) +
                  2 * sizeof(mcpd3::NodeFlow) +
                  3 * sizeof(unsigned char),
          "solver vector estimate should account for arc and node vectors");
  require(estimate.total_bytes ==
              estimate.bk_total_bytes + estimate.solver_vector_bytes,
          "solver total estimate should include BK and solver vectors");
}

mcpd3::PartitionPackage makePackage(const mcpd3::Lagrange &alpha,
                                    const mcpd3::Lagrange &last_alpha) {
  mcpd3::PartitionPackage package;
  package.partition_id = 3;
  package.local_node_count = 2;
  package.arcs = {0, 1};
  package.arc_capacities = {4, 7};
  package.terminal_capacities = {2, -3};
  package.local_to_global = {10, 20};
  package.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/42,
                                        /*global_node_id=*/20,
                                        /*local_index=*/1,
                                        /*is_source=*/true,
                                        /*alpha=*/alpha,
                                        /*last_alpha=*/last_alpha,
                                        /*alpha_momentum=*/0});
  return package;
}

mcpd3::PartitionPackage makeStreamingPackage(int partition_id,
                                             int constraint_id,
                                             int terminal_capacity) {
  mcpd3::PartitionPackage package;
  package.partition_id = partition_id;
  package.local_node_count = 2;
  package.arcs = {0, 1};
  package.arc_capacities = {3 + partition_id, 5 + partition_id};
  package.terminal_capacities = {terminal_capacity, -terminal_capacity - 1};
  package.local_to_global = {100 + partition_id * 10,
                             101 + partition_id * 10};
  package.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/constraint_id,
                                        /*global_node_id=*/101 +
                                            partition_id * 10,
                                        /*local_index=*/1,
                                        /*is_source=*/true,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});
  return package;
}

void requireSolveResultsMatch(const mcpd3::PartitionSolveResult &actual,
                              const mcpd3::PartitionSolveResult &expected,
                              const std::string &context) {
  require(actual.round_id == expected.round_id,
          context + ": round id differs");
  require(actual.partition_id == expected.partition_id,
          context + ": partition id differs");
  require(actual.lower_bound == expected.lower_bound,
          context + ": lower bound differs");
  require(actual.regularization_budget == expected.regularization_budget,
          context + ": regularization budget differs");
  require(actual.regularization_contribution ==
              expected.regularization_contribution,
          context + ": regularization contribution differs");
  require(actual.regularization_anchor_sink_count ==
              expected.regularization_anchor_sink_count,
          context + ": regularization anchor count differs");
  require(actual.regularization_active_sink_count ==
              expected.regularization_active_sink_count,
          context + ": regularization active count differs");
  require(actual.constrained_labels.size() ==
              expected.constrained_labels.size(),
          context + ": constrained label count differs");
  for (size_t i = 0; i < actual.constrained_labels.size(); ++i) {
    require(actual.constrained_labels[i].constraint_id ==
                expected.constrained_labels[i].constraint_id,
            context + ": constraint id differs");
    require(actual.constrained_labels[i].global_node_id ==
                expected.constrained_labels[i].global_node_id,
            context + ": global node id differs");
    require(actual.constrained_labels[i].local_index ==
                expected.constrained_labels[i].local_index,
            context + ": local index differs");
    require(actual.constrained_labels[i].label ==
                expected.constrained_labels[i].label,
            context + ": label differs");
  }
  require(actual.full_labels.size() == expected.full_labels.size(),
          context + ": full label count differs");
  for (size_t i = 0; i < actual.full_labels.size(); ++i) {
    require(actual.full_labels[i].global_node_id ==
                expected.full_labels[i].global_node_id,
            context + ": full label global node differs");
    require(actual.full_labels[i].local_index ==
                expected.full_labels[i].local_index,
            context + ": full label local index differs");
    require(actual.full_labels[i].label == expected.full_labels[i].label,
            context + ": full label value differs");
  }
}

void streamingWorkerMatchesInProcessAcrossEviction() {
  auto package0 = makeStreamingPackage(/*partition_id=*/0,
                                       /*constraint_id=*/10,
                                       /*terminal_capacity=*/2);
  auto package1 = makeStreamingPackage(/*partition_id=*/1,
                                       /*constraint_id=*/11,
                                       /*terminal_capacity=*/4);

  mcpd3::InProcessPartitionWorker reference;
  reference.loadPartition(package0);
  reference.loadPartition(package1);

  mcpd3::StreamingPartitionWorker::Options options;
  options.resident_byte_limit = 1;
  mcpd3::StreamingPartitionWorker streaming(options);
  streaming.loadPartition(package0);
  streaming.loadPartition(package1);

  mcpd3::PartitionSolveRequest first;
  first.round_id = 1;
  first.partition_id = 0;
  first.return_full_labels = true;
  requireSolveResultsMatch(streaming.solveRound(first),
                           reference.solveRound(first),
                           "streaming first solve");

  mcpd3::PartitionSolveRequest second;
  second.round_id = 2;
  second.partition_id = 1;
  second.return_full_labels = true;
  requireSolveResultsMatch(streaming.solveRound(second),
                           reference.solveRound(second),
                           "streaming eviction solve");
  require(streaming.residentPartitionCountForTesting() == 1,
          "streaming worker should keep only one resident partition");
  require(streaming.warmStateWriteCountForTesting() == 1,
          "streaming eviction should persist warm solver state");
  require(streaming.warmStateRestoreCountForTesting() == 0,
          "streaming worker should not restore warm state before reload");

  mcpd3::PartitionSolveRequest third;
  third.round_id = 3;
  third.partition_id = 0;
  third.regularization_strength = 1;
  third.return_full_labels = true;
  third.alpha_updates.push_back(
      mcpd3::AlphaUpdate{/*constraint_id=*/10,
                          /*alpha=*/7,
                          /*last_alpha=*/0,
                          /*alpha_momentum=*/0});
  const auto streaming_third = streaming.solveRound(third);
  const auto reference_third = reference.solveRound(third);
  requireSolveResultsMatch(streaming_third, reference_third,
                           "streaming reload with alpha update");
  require(streaming_third.regularization_budget == 1,
          "first regularized update should consume one budget unit");
  require(streaming.residentPartitionCountForTesting() == 1,
          "streaming reload should preserve the resident budget");
  require(streaming.warmStateWriteCountForTesting() == 2,
          "streaming reload should evict the previous resident warm state");
  require(streaming.warmStateRestoreCountForTesting() == 1,
          "streaming reload should restore persisted warm solver state");

  second.round_id = 4;
  (void)streaming.solveRound(second);
  (void)reference.solveRound(second);

  mcpd3::PartitionSolveRequest fourth = third;
  fourth.round_id = 5;
  fourth.alpha_updates[0].last_alpha = 7;
  fourth.alpha_updates[0].alpha = 8;
  const auto streaming_fourth = streaming.solveRound(fourth);
  const auto reference_fourth = reference.solveRound(fourth);
  requireSolveResultsMatch(streaming_fourth, reference_fourth,
                           "streaming cumulative regularization reload");
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
  require(streaming_fourth.regularization_budget == 1,
          "legacy replay should preserve its binary anchor across reload");
#else
  require(streaming_fourth.regularization_budget == 2,
          "evicted cumulative regularization weight should survive reload");
#endif
}

void streamingWorkerScalesEvictedDiskPayload() {
  auto package0 = makeStreamingPackage(/*partition_id=*/0,
                                       /*constraint_id=*/20,
                                       /*terminal_capacity=*/3);
  auto package1 = makeStreamingPackage(/*partition_id=*/1,
                                       /*constraint_id=*/21,
                                       /*terminal_capacity=*/5);

  mcpd3::InProcessPartitionWorker reference;
  reference.loadPartition(package0);
  reference.loadPartition(package1);

  mcpd3::StreamingPartitionWorker::Options options;
  options.resident_byte_limit = 1;
  mcpd3::StreamingPartitionWorker streaming(options);
  streaming.loadPartition(package0);
  streaming.loadPartition(package1);

  mcpd3::PartitionSolveRequest first;
  first.round_id = 1;
  first.partition_id = 0;
  (void)streaming.solveRound(first);
  (void)reference.solveRound(first);

  mcpd3::PartitionSolveRequest second;
  second.round_id = 2;
  second.partition_id = 1;
  (void)streaming.solveRound(second);
  (void)reference.solveRound(second);
  require(streaming.warmStateWriteCountForTesting() == 1,
          "streaming scale test should persist warm state before scaling");

  streaming.scaleObjective(/*factor=*/2);
  reference.scaleObjective(/*factor=*/2);

  mcpd3::PartitionSolveRequest after_scale;
  after_scale.round_id = 3;
  after_scale.partition_id = 0;
  after_scale.alpha_updates.push_back(
      mcpd3::AlphaUpdate{/*constraint_id=*/20,
                          /*alpha=*/4,
                          /*last_alpha=*/0,
                          /*alpha_momentum=*/0});
  requireSolveResultsMatch(streaming.solveRound(after_scale),
                           reference.solveRound(after_scale),
                           "streaming scaled evicted payload");
  require(streaming.warmStateRestoreCountForTesting() == 0,
          "streaming worker should not restore stale warm state after scaling");
  require(streaming.warmStateWriteCountForTesting() == 2,
          "streaming scaled reload should persist the evicted resident state");
}

void inProcessWorkerScopesObjectiveScalingToSelectedPartitions() {
  auto package0 = makeStreamingPackage(/*partition_id=*/0,
                                       /*constraint_id=*/30,
                                       /*terminal_capacity=*/3);
  auto package1 = makeStreamingPackage(/*partition_id=*/1,
                                       /*constraint_id=*/31,
                                       /*terminal_capacity=*/5);

  mcpd3::InProcessPartitionWorker shared;
  shared.loadPartition(package0);
  shared.loadPartition(package1);

  mcpd3::InProcessPartitionWorker scaled_reference;
  scaled_reference.loadPartition(package0);
  scaled_reference.scaleObjective(/*factor=*/2);
  mcpd3::InProcessPartitionWorker unchanged_reference;
  unchanged_reference.loadPartition(package1);

  bool duplicate_threw = false;
  try {
    shared.scaleObjectivePartitions({0, 0}, /*factor=*/2);
  } catch (const std::runtime_error &) {
    duplicate_threw = true;
  }
  require(duplicate_threw,
          "scoped objective scaling must reject duplicate partition ids");

  bool unknown_threw = false;
  try {
    shared.scaleObjectivePartitions({0, 99}, /*factor=*/2);
  } catch (const std::runtime_error &) {
    unknown_threw = true;
  }
  require(unknown_threw,
          "scoped objective scaling must reject unknown partition ids");

  shared.scaleObjectivePartitions({0}, /*factor=*/2);

  mcpd3::PartitionSolveRequest request0;
  request0.round_id = 1;
  request0.partition_id = 0;
  requireSolveResultsMatch(shared.solveRound(request0),
                           scaled_reference.solveRound(request0),
                           "selected partition objective scaling");

  mcpd3::PartitionSolveRequest request1;
  request1.round_id = 2;
  request1.partition_id = 1;
  requireSolveResultsMatch(shared.solveRound(request1),
                           unchanged_reference.solveRound(request1),
                           "unselected partition objective scaling");

}

void partitionWorkersUnloadSelectedPartitionsAtomically() {
  auto package0 = makeStreamingPackage(/*partition_id=*/0,
                                       /*constraint_id=*/40,
                                       /*terminal_capacity=*/3);
  auto package1 = makeStreamingPackage(/*partition_id=*/1,
                                       /*constraint_id=*/41,
                                       /*terminal_capacity=*/5);

  mcpd3::InProcessPartitionWorker in_process;
  in_process.loadPartition(package0);
  in_process.loadPartition(package1);
  requireThrows([&] { in_process.unloadPartitions({}); },
                "partition unload must reject an empty selection");
  requireThrows([&] { in_process.unloadPartitions({0, 0}); },
                "partition unload must reject duplicate partition ids");
  requireThrows([&] { in_process.unloadPartitions({0, 99}); },
                "partition unload must reject an unknown partition id");

  mcpd3::PartitionSolveRequest request;
  request.partition_id = 0;
  (void)in_process.solveRound(request);
  request.partition_id = 1;
  const auto retained = in_process.solveRound(request);

  in_process.unloadPartitions({0});
  request.partition_id = 0;
  requireThrows([&] { (void)in_process.solveRound(request); },
                "an unloaded in-process partition must become unknown");
  request.partition_id = 1;
  require(in_process.solveRound(request).lower_bound == retained.lower_bound,
          "unloading one partition must preserve the other partition");

  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-unload-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  mcpd3::StreamingPartitionWorker::Options options;
  options.storage_directory = scratch.string();
  options.remove_storage_on_destroy = false;
  mcpd3::StreamingPartitionWorker streaming(options);
  streaming.loadPartition(package0);
  streaming.loadPartition(package1);
  require(std::filesystem::exists(scratch / "partition_0.bin") &&
              std::filesystem::exists(scratch / "partition_1.bin"),
          "streaming unload test requires both disk payloads");
  requireThrows([&] { streaming.unloadPartitions({0, 99}); },
                "streaming unload must validate before deleting payloads");
  require(std::filesystem::exists(scratch / "partition_0.bin") &&
              std::filesystem::exists(scratch / "partition_1.bin"),
          "invalid streaming unload must not partially delete payloads");

  streaming.unloadPartitions({0});
  require(!std::filesystem::exists(scratch / "partition_0.bin") &&
              std::filesystem::exists(scratch / "partition_1.bin"),
          "streaming unload must remove only the selected disk payload");
  request.partition_id = 0;
  requireThrows([&] { (void)streaming.solveRound(request); },
                "an unloaded streaming partition must become unknown");
  request.partition_id = 1;
  (void)streaming.solveRound(request);
  streaming.unloadPartitions({1});
  require(!std::filesystem::exists(scratch / "partition_1.bin"),
          "streaming unload must remove the final disk payload");
  std::filesystem::remove_all(scratch);
}

struct DirectSolverResult {
  mcpd3::Objective lower_bound;
  int constrained_label;
  mcpd3::Objective regularization_budget;
  mcpd3::Objective regularization_contribution;
  long regularization_anchor_sink_count;
  long regularization_active_sink_count;
};

DirectSolverResult solveDirect(mcpd3::DualDecompositionConstraintArcReference ref,
                               mcpd3::PrimalDualMinCutSolver *solver,
                               const mcpd3::Capacity &regularization_strength) {
  solver->setRegularizationStrength(regularization_strength);
  solver->solve();
  return DirectSolverResult{solver->getMinCutValue(),
                            solver->getMinCutSolution(ref->local_index_source),
                            solver->getLastRegularizationBudget(),
                            solver->getLastRegularizationContribution(),
                            solver->getLastRegularizationAnchorSinkCount(),
                            solver->getLastRegularizationActiveSinkCount()};
}

void requireMatchesDirect(const mcpd3::PartitionSolveResult &worker_result,
                          const DirectSolverResult &direct_result,
                          long round_id) {
  require(worker_result.round_id == round_id, "round id was not preserved");
  require(worker_result.partition_id == 3, "partition id was not preserved");
  require(worker_result.lower_bound == direct_result.lower_bound,
          "worker lower bound differs from direct solver");
  require(worker_result.regularization_budget ==
              direct_result.regularization_budget,
          "regularization budget differs from direct solver");
  require(worker_result.regularization_contribution ==
              direct_result.regularization_contribution,
          "regularization contribution differs from direct solver");
  require(worker_result.regularization_anchor_sink_count ==
              direct_result.regularization_anchor_sink_count,
          "regularization anchor sink count differs from direct solver");
  require(worker_result.regularization_active_sink_count ==
              direct_result.regularization_active_sink_count,
          "regularization active sink count differs from direct solver");
  require(worker_result.constrained_labels.size() == 1,
          "expected one constrained label");
  require(worker_result.constrained_labels[0].constraint_id == 42,
          "constraint id was not preserved");
  require(worker_result.constrained_labels[0].global_node_id == 20,
          "global node id was not preserved");
  require(worker_result.constrained_labels[0].local_index == 1,
          "local index was not preserved");
  require(worker_result.constrained_labels[0].label ==
              direct_result.constrained_label,
          "worker constrained label differs from direct solver");
}

void inProcessPartitionWorkerMatchesDirectSolverAcrossAlphaUpdate() {
  auto package = makePackage(/*alpha=*/0, /*last_alpha=*/0);

  std::list<mcpd3::DualDecompositionConstraintArc> direct_constraints;
  direct_constraints.emplace_back(/*alpha=*/0, /*last_alpha=*/0,
                                  /*alpha_momentum=*/0,
                                  /*partition_index_source=*/3,
                                  /*partition_index_target=*/4,
                                  /*local_index_source=*/1,
                                  /*local_index_target=*/-1);
  auto direct_ref = --direct_constraints.end();
  mcpd3::PrimalDualMinCutSolver direct_solver(
      package.local_node_count,
      static_cast<int>(package.arcs.size() / 2), package.arcs.copyToVector(),
      package.arc_capacities.copyToVector(),
      package.terminal_capacities.copyToVector());
  direct_solver.addSourceDualDecompositionConstraint(direct_ref);

  mcpd3::InProcessPartitionWorker worker;
  worker.loadPartition(package);

  mcpd3::PartitionSolveRequest first_request;
  first_request.round_id = 100;
  first_request.scale = 10000;
  first_request.regularization_strength = 0;
  auto direct_first = solveDirect(direct_ref, &direct_solver,
                                  first_request.regularization_strength);
  auto worker_first = worker.solveRound(first_request);
  requireMatchesDirect(worker_first, direct_first, first_request.round_id);

  direct_ref->last_alpha = direct_ref->alpha;
  direct_ref->alpha = 6;
  mcpd3::PartitionSolveRequest second_request;
  second_request.round_id = 101;
  second_request.scale = 10000;
  second_request.regularization_strength = 1;
  second_request.alpha_updates.push_back(
      mcpd3::AlphaUpdate{/*constraint_id=*/42,
                          /*alpha=*/direct_ref->alpha,
                          /*last_alpha=*/direct_ref->last_alpha,
                          /*alpha_momentum=*/direct_ref->alpha_momentum});

  auto direct_second = solveDirect(direct_ref, &direct_solver,
                                   second_request.regularization_strength);
  auto worker_second = worker.solveRound(second_request);
  requireMatchesDirect(worker_second, direct_second, second_request.round_id);
}

void inProcessPartitionWorkerReturnsFullLabelsOnRequest() {
  auto package = makePackage(/*alpha=*/0, /*last_alpha=*/0);

  mcpd3::InProcessPartitionWorker worker;
  worker.loadPartition(package);

  mcpd3::PartitionSolveRequest default_request;
  default_request.round_id = 1;
  default_request.partition_id = package.partition_id;
  const auto default_result = worker.solveRound(default_request);
  require(default_result.full_labels.empty(),
          "worker should not return full labels unless requested");

  mcpd3::PartitionSolveRequest full_request;
  full_request.round_id = 2;
  full_request.partition_id = package.partition_id;
  full_request.return_full_labels = true;
  const auto full_result = worker.solveRound(full_request);
  require(full_result.full_labels.size() == 2,
          "worker should return one full label per local node");
  require(full_result.full_labels[0].global_node_id == 10,
          "first full label global node id differs");
  require(full_result.full_labels[0].local_index == 0,
          "first full label local index differs");
  require(full_result.full_labels[1].global_node_id == 20,
          "second full label global node id differs");
  require(full_result.full_labels[1].local_index == 1,
          "second full label local index differs");
  require(full_result.full_labels[0].label == worker.minCutSolution(3)[0] &&
              full_result.full_labels[1].label == worker.minCutSolution(3)[1],
          "full labels should match the worker min-cut solution");

  std::vector<mcpd3::NodeLabel> copied_labels(2);
  worker.copyFullLabels(package.partition_id, /*offset=*/0,
                        copied_labels.data(),
                        copied_labels.size());
  require(copied_labels[0].global_node_id == 10 &&
              copied_labels[0].local_index == 0 &&
              copied_labels[0].label == full_result.full_labels[0].label &&
              copied_labels[1].global_node_id == 20 &&
              copied_labels[1].local_index == 1 &&
              copied_labels[1].label == full_result.full_labels[1].label,
          "bounded label recovery should copy the current solved labeling");
  requireThrows(
      [&] {
        worker.copyFullLabels(package.partition_id, /*offset=*/2,
                              copied_labels.data(), 1);
      },
      "bounded label recovery should reject an out-of-range destination");
  requireThrows(
      [&] {
        worker.copyFullLabels(package.partition_id, /*offset=*/0, nullptr, 2);
      },
      "bounded label recovery should reject a null destination");
  requireThrows(
      [&] {
        worker.copyFullLabels(/*partition_id=*/999, /*offset=*/0, nullptr, 0);
      },
      "bounded label recovery should reject an unknown partition");
  require(worker.fullLabelCount(package.partition_id) == 2,
          "worker should report its exact full label count");
  mcpd3::NodeLabel second_label;
  worker.copyFullLabels(package.partition_id, /*offset=*/1, &second_label, 1);
  require(second_label.global_node_id == 20 && second_label.local_index == 1 &&
              second_label.label == full_result.full_labels[1].label,
          "bounded label recovery should support partial ranges");
}

long countWorkerDisagreements(
    const std::vector<mcpd3::PartitionSolveResult> &results) {
  std::map<int, std::vector<int>> labels_by_constraint;
  for (const auto &result : results) {
    for (const auto &label : result.constrained_labels) {
      labels_by_constraint[label.constraint_id].push_back(label.label);
    }
  }

  long disagreement_count = 0;
  for (const auto &entry : labels_by_constraint) {
    const auto &labels = entry.second;
    require(labels.size() == 2, "expected pairwise constraint labels");
    disagreement_count += std::abs(labels[1] - labels[0]);
  }
  return disagreement_count;
}

void exportedPartitionPackagesMatchDualDecompositionRound() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.canonical_cut_selection =
      mcpd3::CanonicalCutSelection::MAXIMUM_LABELS;
  options.force_full_mincut_recompute = true;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  const auto &packages = dual_decomp.getPartitionPackages();
  require(packages.size() == 2, "expected two exported partition packages");

  std::vector<mcpd3::PartitionSolveResult> worker_results;
  mcpd3::Objective worker_lower_bound = 0;
  for (const auto &package : packages) {
    mcpd3::InProcessPartitionWorker worker;
    worker.loadPartition(package);

    mcpd3::PartitionSolveRequest request;
    request.round_id = 1;
    request.scale = 100;
    request.regularization_strength = 0;
    auto result = worker.solveRound(request);
    worker_lower_bound += result.lower_bound;
    worker_results.push_back(result);
  }

  dual_decomp.runOptimizationScale(
      /*nstep=*/1, /*step_size=*/100, /*max_cycle_count=*/2,
      /*use_momentum=*/false);

  require(worker_lower_bound == dual_decomp.getBestLowerBoundRaw(),
          "worker lower bound differs from DualDecomposition");
  require(countWorkerDisagreements(worker_results) ==
              dual_decomp.getLastDisagreementCount(),
          "worker disagreement count differs from DualDecomposition");
}

void nativeProgressReportsEveryPartitionSolve() {
  const char *old_progress = std::getenv("MCPD3_PROGRESS");
  const bool had_old_progress = old_progress != nullptr;
  const std::string old_progress_value = had_old_progress ? old_progress : "";
  ::unsetenv("MCPD3_PROGRESS");
  ::setenv("MCPD3_PARTITIONER", "basic", 1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 100;
  options.force_full_mincut_recompute = true;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  ::setenv("MCPD3_PROGRESS", "1", 1);
  const std::string progress = captureStderr([&] {
    dual_decomp.runOptimizationScale(
        /*nstep=*/1, /*step_size=*/100, /*max_cycle_count=*/2,
        /*use_momentum=*/false);
  });

  if (had_old_progress) {
    ::setenv("MCPD3_PROGRESS", old_progress_value.c_str(), 1);
  } else {
    ::unsetenv("MCPD3_PROGRESS");
  }

  const size_t start_zero = progress.find(
      "stage=dd_partition_solve_started objective_scale=1 iter=0 partition=0 ");
  const size_t complete_zero = progress.find(
      "stage=dd_partition_solve_completed objective_scale=1 iter=0 partition=0 ");
  const size_t start_one = progress.find(
      "stage=dd_partition_solve_started objective_scale=1 iter=0 partition=1 ");
  const size_t complete_one = progress.find(
      "stage=dd_partition_solve_completed objective_scale=1 iter=0 partition=1 ");
  if (start_zero == std::string::npos ||
      complete_zero == std::string::npos || start_one == std::string::npos ||
      complete_one == std::string::npos) {
    std::cerr << "captured partition progress:\n" << progress;
  }
  require(start_zero != std::string::npos &&
              complete_zero != std::string::npos &&
              start_one != std::string::npos &&
              complete_one != std::string::npos,
          "native progress must report each partition start and completion");
  require(start_zero < complete_zero && complete_zero < start_one &&
              start_one < complete_one,
          "a one-thread solve must expose sequential partition handoff");
  require(progress.find("partition_count=2") != std::string::npos,
          "partition progress must report the round partition count");
  require(progress.find("solve_elapsed_sec=") != std::string::npos,
          "partition completion must report local maxflow wall time");
}

void exportedPartitionPackagesMaterializeBoundaryDuplicates() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.emit_partition_packages = true;
  options.construct_solvers = false;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/5,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{1, 3},
      /*arc_capacities=*/std::vector<int>{6, 6},
      /*terminal_capacities=*/std::vector<int>{0, 0, 0, 7, 0}, options);

  const auto &packages = dual_decomp.getPartitionPackages();
  require(packages.size() == 2, "expected two exported partition packages");
  require(packages[0].local_node_count == 2,
          "owning partition should contain both arc endpoints");
  require(packages[0].local_to_global == std::vector<int>({1, 3}),
          "owning partition local_to_global should preserve arc endpoints");
  require(packages[0].terminal_capacities ==
              std::vector<mcpd3::Capacity>({0, 0}),
          "arc-owner boundary clone must not receive the terminal");
  require(packages[1].local_node_count == 1,
          "target partition should contain an isolated boundary duplicate");
  require(packages[1].terminal_capacities.size() == 1,
          "target duplicate should have a zero terminal entry");
  require(packages[1].terminal_capacities[0] == 7,
          "boundary node home partition should receive the terminal once");
  require(packages[1].local_to_global == std::vector<int>({3}),
          "target duplicate should map back to the constrained global node");
  require(packages[0].constraint_endpoints.size() == 1 &&
              packages[1].constraint_endpoints.size() == 1,
          "cross-partition arc should emit one constraint endpoint per side");
}

void disabledPartitionPackageExportPreservesNativeSolve() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.emit_partition_packages = false;
  options.use_momentum = false;
  options.enable_group_stopping = false;

  mcpd3::DualDecomposition native_only(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  bool package_access_threw = false;
  try {
    (void)native_only.getPartitionPackages();
  } catch (const std::runtime_error &e) {
    package_access_threw =
        std::string(e.what()).find("partition package export is disabled") !=
        std::string::npos;
  }
  require(package_access_threw,
          "disabled partition package export should reject package access");

  options.emit_partition_packages = true;
  mcpd3::DualDecomposition exported(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  native_only.runOptimizationScale(
      /*nstep=*/1, /*step_size=*/100, /*max_cycle_count=*/2,
      /*use_momentum=*/false);
  exported.runOptimizationScale(
      /*nstep=*/1, /*step_size=*/100, /*max_cycle_count=*/2,
      /*use_momentum=*/false);

  require(native_only.getBestLowerBoundRaw() == exported.getBestLowerBoundRaw(),
          "native-only lower bound should match package-export solve");
  require(native_only.getLastDisagreementCount() ==
              exported.getLastDisagreementCount(),
          "native-only disagreement count should match package-export solve");
  require(native_only.getTotalOptimizationIterations() == 1,
          "native-only solve should run without exported packages");
}

void requirePackagesEqual(const mcpd3::PartitionPackage &lhs,
                          const mcpd3::PartitionPackage &rhs) {
  require(lhs.partition_id == rhs.partition_id, "partition id differs");
  require(lhs.local_node_count == rhs.local_node_count,
          "local node count differs");
  require(lhs.arcs == rhs.arcs, "package arcs differ");
  require(lhs.arc_capacities == rhs.arc_capacities,
          "package arc capacities differ");
  require(lhs.terminal_capacities == rhs.terminal_capacities,
          "package terminal capacities differ");
  require(lhs.local_to_global == rhs.local_to_global,
          "package local_to_global differs");
  require(lhs.objective_multiplier == rhs.objective_multiplier,
          "package objective multiplier differs");
  require(lhs.canonical_cut_selection == rhs.canonical_cut_selection,
          "package canonical cut selection differs");
  require(lhs.force_full_mincut_recompute == rhs.force_full_mincut_recompute,
          "package full recompute policy differs");
  require(lhs.reference_cut_labels == rhs.reference_cut_labels,
          "package reference cut labels differ");
  require(lhs.reference_cut_selection == rhs.reference_cut_selection,
          "package reference cut selection differs");
  require(lhs.reference_cut_check_interval ==
              rhs.reference_cut_check_interval,
          "package reference cut interval differs");
  require(lhs.constraint_endpoints.size() ==
              rhs.constraint_endpoints.size(),
          "constraint endpoint count differs");
  for (size_t i = 0; i < lhs.constraint_endpoints.size(); ++i) {
    const auto &left = lhs.constraint_endpoints[i];
    const auto &right = rhs.constraint_endpoints[i];
    require(left.constraint_id == right.constraint_id,
            "constraint endpoint id differs");
    require(left.global_node_id == right.global_node_id,
            "constraint endpoint global node differs");
    require(left.local_index == right.local_index,
            "constraint endpoint local index differs");
    require(left.is_source == right.is_source,
            "constraint endpoint side differs");
    require(left.alpha == right.alpha, "constraint endpoint alpha differs");
    require(left.last_alpha == right.last_alpha,
            "constraint endpoint last alpha differs");
    require(left.alpha_momentum == right.alpha_momentum,
            "constraint endpoint alpha momentum differs");
  }
}

void haloDepthOneExportsIdenticalLegacyPackages() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.emit_partition_packages = true;
  options.construct_solvers = false;

  const std::vector<int> arcs{0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
  const std::vector<int> capacities{2, 3, 4, 5, 6,
                                    7, 8, 9, 10, 11};
  const std::vector<int> terminals{0, 2, 0, -3, 0, 0};
  mcpd3::DualDecomposition legacy(2, 6, 5, arcs, capacities, terminals,
                                  options);

  auto halo_options = options;
  halo_options.halo_depth = 1;
  mcpd3::DualDecomposition halo(2, 6, 5, arcs, capacities, terminals,
                                halo_options);

  const auto &legacy_packages = legacy.getPartitionPackages();
  const auto &halo_packages = halo.getPartitionPackages();
  require(halo.getHaloObjectiveMultiplier() == 1,
          "h1 must retain a unit halo objective multiplier");
  require(legacy_packages.size() == halo_packages.size(),
          "h1 package count must match the legacy decomposition");
  for (size_t partition = 0; partition < legacy_packages.size(); ++partition) {
    requirePackagesEqual(legacy_packages[partition],
                         halo_packages[partition]);
  }
}

void haloDepthTwoExportsScaledInducedSubproblems() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.emit_partition_packages = true;
  options.construct_solvers = false;
  options.halo_depth = 2;

  mcpd3::DualDecomposition halo(
      /*npartition=*/2, /*nnode=*/6, /*narc=*/5,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3, 3, 4, 4, 5},
      /*arc_capacities=*/std::vector<int>{1, 3, 1, 3, 1, 3, 1, 3, 1, 3},
      /*terminal_capacities=*/std::vector<int>{6, 4, 0, 0, 0, -8}, options);

  require(halo.getHaloObjectiveMultiplier() == 2,
          "h2 path decomposition should expose multiplier two");
  const auto &packages = halo.getPartitionPackages();
  require(packages.size() == 2, "h2 path should export two packages");
  require(packages[0].objective_multiplier == 2 &&
              packages[1].objective_multiplier == 2,
          "every h2 package should carry the common objective multiplier");
  require(packages[0].local_to_global ==
              std::vector<int>({0, 1, 2, 3, 4}),
          "left h2 package should contain its two-step node halo");
  require(packages[0].arcs ==
              std::vector<int>({0, 1, 1, 2, 2, 3, 3, 4}),
          "left h2 package should contain every induced halo edge");
  require(packages[0].arc_capacities ==
              std::vector<mcpd3::Capacity>({2, 6, 1, 3, 1, 3, 1, 3}),
          "left h2 edge capacities should use Q divided by multiplicity");
  require(packages[0].terminal_capacities ==
              std::vector<mcpd3::Capacity>({12, 4, 0, 0, 0}),
          "left h2 unaries should use Q divided by node multiplicity");

  require(packages[1].local_to_global ==
              std::vector<int>({1, 2, 3, 4, 5}),
          "right h2 package should contain its two-step node halo");
  require(packages[1].arcs ==
              std::vector<int>({0, 1, 1, 2, 2, 3, 3, 4}),
          "right h2 package should contain every induced halo edge");
  require(packages[1].arc_capacities ==
              std::vector<mcpd3::Capacity>({1, 3, 1, 3, 1, 3, 2, 6}),
          "right h2 edge capacities should use Q divided by multiplicity");
  require(packages[1].terminal_capacities ==
              std::vector<mcpd3::Capacity>({4, 0, 0, 0, -16}),
          "right h2 unaries should use Q divided by node multiplicity");

  for (const auto &package : packages) {
    require(package.constraint_endpoints.size() == 4,
            "every duplicated h2 node should have a consensus endpoint");
    std::vector<int> constrained_nodes;
    for (const auto &endpoint : package.constraint_endpoints) {
      constrained_nodes.push_back(endpoint.global_node_id);
    }
    std::sort(constrained_nodes.begin(), constrained_nodes.end());
    require(constrained_nodes == std::vector<int>({1, 2, 3, 4}),
            "h2 constraints should cover all and only duplicated nodes");
  }
}

void haloDepthTwoCertifiesTheUnscaledDirectObjective() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  const std::vector<int> arcs{0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
  const std::vector<int> capacities{3, 3, 5, 5, 2,
                                    2, 7, 7, 4, 4};
  const std::vector<int> terminals{20, 0, 0, 0, 0, -20};

  mcpd3::PrimalDualMinCutSolver direct(
      /*nnode=*/6, /*narc=*/5, std::vector<int>(arcs), capacities, terminals);
  direct.solve();
  const mcpd3::Objective exact = direct.getMinCutValue();

  mcpd3::DualDecompositionOptions options;
  options.verbose = false;
  options.thread_count = 1;
  options.halo_depth = 2;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 200;
  options.patience = 200;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;

  mcpd3::DualDecomposition halo(2, 6, 5, arcs, capacities, terminals,
                                options);
  halo.solve();

  require(halo.getLastDisagreementCount() == 0,
          "h2 path solve should reach primal agreement");
  require(halo.getBestCertifiedLowerBoundRaw() ==
              exact * halo.getHaloObjectiveMultiplier(),
          "h2 raw certificate should equal Q times the direct objective");
  require(halo.getBestCertifiedLowerBound() ==
              mcpd3::integer_to_double(exact),
          "h2 public certificate should normalize away the halo multiplier");
}

void haloDepthTwoPromotesWhenDuplicatedNodesExhaustTheBudget() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  const std::vector<int> arcs{0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
  const std::vector<int> capacities{10, 10, 10, 10, 10,
                                    10, 10, 10, 10, 10};
  const std::vector<int> terminals{100, 0, 0, 0, 0, -100};

  mcpd3::PrimalDualMinCutSolver direct(
      /*nnode=*/6, /*narc=*/5, std::vector<int>(arcs), capacities, terminals);
  direct.solve();

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.halo_depth = 2;
  options.objective_scale = 10;
  options.initial_step_size = 10;
  options.num_optimization_scales = 3;
  options.max_iteration_count = 100;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.max_objective_scale_promotions = 2;

  mcpd3::DualDecomposition halo(
      /*npartition=*/2, /*nnode=*/6, /*narc=*/5,
      arcs, capacities, terminals, options);
  require(halo.getHaloObjectiveMultiplier() == 2,
          "h2 promotion fixture must exercise a nonunit halo multiplier");

  halo.solve();

  require(halo.getObjectiveScalePromotionCount() == 1,
          "h2 duplicated-node regularization should force one promotion");
  require(halo.getScale() == 200,
          "h2 promotion should retain Q=2 while increasing the primary "
          "scale by one decade");
  require(halo.getLastRegularizationBudget() < halo.getScale(),
          "promoted h2 solve should finish under the strict scaled budget");
  require(halo.getLastDisagreementCount() == 0,
          "promoted h2 solve should preserve state and reach agreement");
  require(halo.getBestCertifiedLowerBoundRaw() ==
              direct.getMinCutValue() * halo.getHaloObjectiveMultiplier() *
                  10,
          "promoted h2 agreement must certify the exactly scaled direct "
          "objective");
  require(halo.getBestCertifiedLowerBound() ==
              mcpd3::integer_to_double(direct.getMinCutValue()) / 10.0,
          "public h2 objective must normalize both Q and primary scaling");
}

template <typename ArcContainer, typename CapacityContainer,
          typename TerminalContainer>
mcpd3::Objective evaluateCutObjective(
    const ArcContainer &arcs, const CapacityContainer &capacities,
    const TerminalContainer &terminals,
    const std::vector<int> &labels) {
  mcpd3::Objective objective = 0;
  for (size_t arc = 0; arc < arcs.size() / 2; ++arc) {
    const int source = arcs[2 * arc];
    const int target = arcs[2 * arc + 1];
    if (labels[source] == 0 && labels[target] == 1) {
      objective = mcpd3::checked_add(
          objective, mcpd3::widen_capacity(capacities[2 * arc]));
    } else if (labels[source] == 1 && labels[target] == 0) {
      objective = mcpd3::checked_add(
          objective, mcpd3::widen_capacity(capacities[2 * arc + 1]));
    }
  }
  for (size_t node = 0; node < terminals.size(); ++node) {
    if (labels[node] == 0 && terminals[node] < 0) {
      objective = mcpd3::checked_add(
          objective, mcpd3::absolute_capacity(terminals[node]));
    } else if (labels[node] == 1 && terminals[node] > 0) {
      objective = mcpd3::checked_add(
          objective, mcpd3::widen_capacity(terminals[node]));
    }
  }
  return objective;
}

void haloLocalObjectivesSumToScaledGlobalObjectiveExhaustively() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  std::mt19937 generator(1873);
  std::uniform_int_distribution<int> capacity_distribution(0, 7);
  std::uniform_int_distribution<int> terminal_distribution(-5, 5);
  const std::vector<int> depths{1, 2, 3, mcpd3::kInfiniteHaloDepth};

  for (int graph_index = 0; graph_index < 8; ++graph_index) {
    const int node_count = 6;
    std::vector<int> arcs{0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
    if (graph_index % 2 == 0) {
      arcs.insert(arcs.end(), {0, 3, 1, 4, 2, 5});
    }
    std::vector<mcpd3::Capacity> capacities(arcs.size());
    for (auto &capacity : capacities) {
      capacity = capacity_distribution(generator);
    }
    std::vector<mcpd3::Capacity> terminals(node_count);
    for (auto &terminal : terminals) {
      terminal = terminal_distribution(generator);
    }

    for (const int depth : depths) {
      mcpd3::DualDecompositionOptions options;
      options.track_primal_upper_bound = false;
      options.verbose = false;
      options.thread_count = 1;
      options.emit_partition_packages = true;
      options.construct_solvers = false;
      options.halo_depth = depth;
      mcpd3::DualDecomposition halo(
          /*npartition=*/3, node_count,
          static_cast<int>(arcs.size() / 2), arcs, capacities, terminals,
          options);

      for (int mask = 0; mask < (1 << node_count); ++mask) {
        std::vector<int> global_labels(static_cast<size_t>(node_count));
        for (int node = 0; node < node_count; ++node) {
          global_labels[static_cast<size_t>(node)] = (mask >> node) & 1;
        }
        const mcpd3::Objective expected = mcpd3::checked_scale(
            evaluateCutObjective(arcs, capacities, terminals, global_labels),
            halo.getHaloObjectiveMultiplier());
        mcpd3::Objective actual = 0;
        for (const auto &package : halo.getPartitionPackages()) {
          std::vector<int> local_labels(package.local_to_global.size());
          for (size_t local = 0; local < local_labels.size(); ++local) {
            local_labels[local] = global_labels[static_cast<size_t>(
                package.local_to_global[local])];
          }
          actual = mcpd3::checked_add(
              actual,
              evaluateCutObjective(package.arcs, package.arc_capacities,
                                   package.terminal_capacities, local_labels));
        }
        require(actual == expected,
                "agreed halo local objectives must sum to Q times the global "
                "objective");
      }
    }
  }
}

void haloCapacityReplacementPreservesMultiplicityScaling() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  const std::vector<int> arcs{0, 1, 1, 2, 2, 3, 3, 4, 4, 5};
  const std::vector<mcpd3::Capacity> initial_capacities{
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  const std::vector<mcpd3::Capacity> initial_terminals{10, 0, 0, 0, 0, -10};
  const std::vector<mcpd3::Capacity> replacement_capacities{
      2, 4, 3, 5, 4, 6, 5, 7, 6, 8};
  const std::vector<mcpd3::Capacity> replacement_terminals{
      9, 3, 0, 0, -2, -11};

  mcpd3::DualDecompositionOptions options;
  options.verbose = false;
  options.thread_count = 1;
  options.halo_depth = 2;
  options.emit_partition_packages = true;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 300;
  options.patience = 300;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;

  mcpd3::DualDecomposition halo(2, 6, 5, arcs, initial_capacities,
                                initial_terminals, options);
  halo.replaceProblemCapacities(
      replacement_capacities, replacement_terminals,
      /*preserve_alpha_state=*/false, /*preserve_flow_state=*/false);

  const auto &packages = halo.getPartitionPackages();
  require(packages[0].arc_capacities ==
              std::vector<mcpd3::Capacity>({4, 8, 3, 5, 4, 6, 5, 7}),
          "left replacement arcs must preserve halo multiplicity factors");
  require(packages[1].arc_capacities ==
              std::vector<mcpd3::Capacity>({3, 5, 4, 6, 5, 7, 12, 16}),
          "right replacement arcs must preserve halo multiplicity factors");
  require(packages[0].terminal_capacities ==
              std::vector<mcpd3::Capacity>({18, 3, 0, 0, -2}),
          "left replacement unaries must preserve halo multiplicity factors");
  require(packages[1].terminal_capacities ==
              std::vector<mcpd3::Capacity>({3, 0, 0, -2, -22}),
          "right replacement unaries must preserve halo multiplicity factors");

  mcpd3::PrimalDualMinCutSolver direct(
      6, 5, std::vector<int>(arcs), replacement_capacities,
      replacement_terminals);
  direct.solve();
  halo.scaleProblem(3);
  halo.solve();
  require(halo.getLastDisagreementCount() == 0,
          "replacement h2 solve should reach agreement");
  require(halo.getBestCertifiedLowerBoundRaw() ==
              direct.getMinCutValue() * halo.getHaloObjectiveMultiplier() * 3,
          "promoted replacement h2 certificate should retain exact Q scaling");
  require(halo.getBestCertifiedLowerBound() ==
              mcpd3::integer_to_double(direct.getMinCutValue()),
          "promotion should preserve the normalized replacement objective");
}

void haloPackageCoordinatorNormalizesTheObjectiveMultiplier() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  const std::vector<int> arcs{0, 1, 1, 2, 2, 3};
  const std::vector<int> capacities{3, 3, 4, 4, 5, 5};
  const std::vector<int> terminals{20, 0, 0, -20};

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.halo_depth = mcpd3::kInfiniteHaloDepth;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 10;
  options.patience = 10;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;

  mcpd3::DualDecomposition native(2, 4, 3, arcs, capacities, terminals,
                                  options);
  native.solve();

  auto export_options = options;
  export_options.construct_solvers = false;
  mcpd3::DualDecomposition package_source(
      2, 4, 3, arcs, capacities, terminals, export_options);
  std::vector<mcpd3::PartitionPackage> packages =
      package_source.getPartitionPackages();
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());

  mcpd3::PartitionWorkerCoordinatorOptions coordinator_options;
  coordinator_options.num_optimization_scales = 1;
  coordinator_options.max_iteration_count = 10;
  coordinator_options.initial_step_size = 1;
  coordinator_options.patience = 10;
  coordinator_options.use_momentum = false;
  coordinator_options.enable_group_stopping = false;
  coordinator_options.regularization_scheme =
      mcpd3::PartitionWorkerRegularizationScheme::NONE;
  mcpd3::PartitionWorkerCoordinator coordinator(
      std::move(packages), std::move(workers), coordinator_options);
  const auto result = coordinator.solve();

  require(result.scale == 2,
          "package coordinator scale should include the halo multiplier");
  require(result.final_disagreement_count == 0,
          "identical infinite-halo packages should agree immediately");
  require(result.best_certified_lower_bound_raw ==
              native.getBestCertifiedLowerBoundRaw(),
          "package and native halo raw certificates should match");
  require(result.best_certified_lower_bound ==
              native.getBestCertifiedLowerBound(),
          "package and native halo normalized certificates should match");

  auto mismatched_packages = package_source.getPartitionPackages();
  mismatched_packages[1].objective_multiplier = 3;
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> invalid_workers;
  invalid_workers.push_back(
      std::make_unique<mcpd3::InProcessPartitionWorker>());
  invalid_workers.push_back(
      std::make_unique<mcpd3::InProcessPartitionWorker>());
  requireThrows(
      [&] {
        mcpd3::PartitionWorkerCoordinator invalid(
            std::move(mismatched_packages), std::move(invalid_workers),
            coordinator_options);
      },
      "coordinator should reject mismatched halo objective multipliers");

  auto invalid_package = package_source.getPartitionPackages().front();
  invalid_package.objective_multiplier = 0;
  requireThrows(
      [&] { mcpd3::validatePartitionPackage(invalid_package); },
      "worker package validation should reject a zero objective multiplier");
}

#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
void haloFlowHeatAggregatesEveryLocalEdgeCopy() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.halo_depth = mcpd3::kInfiniteHaloDepth;
  options.track_arc_flow_updates = true;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 1;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;

  mcpd3::DualDecomposition halo(
      /*npartition=*/2, /*nnode=*/4, /*narc=*/3,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3},
      /*arc_capacities=*/std::vector<int>{5, 5, 5, 5, 5, 5},
      /*terminal_capacities=*/std::vector<int>{20, 0, 0, -20}, options);
  halo.runOptimizationScale(1, 1, 2, false);

  const auto counts = halo.getArcFlowUpdateCounts();
  require(counts == std::vector<std::uint64_t>({2, 2, 2}),
          "global halo flow heat should sum both full-graph edge copies");
  halo.resetArcFlowUpdateCounts();
  require(halo.getArcFlowUpdateCounts() ==
              std::vector<std::uint64_t>({0, 0, 0}),
          "halo flow-heat reset should clear every local edge copy");
}
#endif

void packageOnlyExportMatchesSolverBackedExport() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.emit_partition_packages = true;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.force_full_mincut_recompute = true;
  options.reference_cut_labels = {1, 0};
  options.reference_cut_selection =
      mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL;
  options.reference_cut_check_interval = 2;

  mcpd3::DualDecomposition solver_backed(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  auto package_only_options = options;
  package_only_options.construct_solvers = false;
  mcpd3::DualDecomposition package_only(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, package_only_options);

  const auto &solver_packages = solver_backed.getPartitionPackages();
  const auto &package_only_packages = package_only.getPartitionPackages();
  require(solver_packages.size() == package_only_packages.size(),
          "package-only export count should match solver-backed export");
  for (size_t i = 0; i < solver_packages.size(); ++i) {
    requirePackagesEqual(solver_packages[i], package_only_packages[i]);
  }

  bool solve_threw = false;
  try {
    package_only.solve();
  } catch (const std::runtime_error &e) {
    solve_threw =
        std::string(e.what()).find("requires constructed solvers") !=
        std::string::npos;
  }
  require(solve_threw, "package-only DualDecomposition should reject solve");

  auto invalid_options = package_only_options;
  invalid_options.emit_partition_packages = false;
  bool invalid_threw = false;
  try {
    mcpd3::DualDecomposition invalid(
        /*npartition=*/2,
        /*nnode=*/2,
        /*narc=*/1,
        /*arcs=*/std::vector<int>{0, 1},
        /*arc_capacities=*/std::vector<int>{3, 5},
        /*terminal_capacities=*/std::vector<int>{2, -4}, invalid_options);
  } catch (const std::runtime_error &e) {
    invalid_threw =
        std::string(e.what()).find("solver construction") !=
        std::string::npos;
  }
  require(invalid_threw,
          "package-only construction should require package export");
}

mcpd3::DualDecomposition makeTinyDualDecomposition() {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;

  return mcpd3::DualDecomposition(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
}

void requireConstraintSnapshotsEqual(
    const std::vector<mcpd3::DualDecompositionConstraintSnapshot> &actual,
    const std::vector<mcpd3::DualDecompositionConstraintSnapshot> &expected,
    const std::string &context) {
  require(actual.size() == expected.size(),
          context + ": alpha snapshot count differs");
  for (size_t i = 0; i < actual.size(); ++i) {
    const auto &lhs = actual[i];
    const auto &rhs = expected[i];
    const std::string item_context =
        context + ": constraint " + std::to_string(i);
    require(lhs.constraint_id == rhs.constraint_id,
            item_context + " id differs");
    require(lhs.global_node_id == rhs.global_node_id,
            item_context + " global node differs");
    require(lhs.partition_index_source == rhs.partition_index_source,
            item_context + " source partition differs");
    require(lhs.partition_index_target == rhs.partition_index_target,
            item_context + " target partition differs");
    require(lhs.local_index_source == rhs.local_index_source,
            item_context + " source local index differs");
    require(lhs.local_index_target == rhs.local_index_target,
            item_context + " target local index differs");
    require(lhs.alpha == rhs.alpha, item_context + " alpha differs");
    require(lhs.last_alpha == rhs.last_alpha,
            item_context + " last alpha differs");
    require(lhs.alpha_momentum == rhs.alpha_momentum,
            item_context + " alpha momentum differs");
  }
}

std::vector<mcpd3::DualDecompositionPartitionSnapshot>
partitionSnapshotsFromTrace(const mcpd3::PartitionWorkerRoundTrace &trace) {
  std::vector<mcpd3::DualDecompositionPartitionSnapshot> snapshots;
  snapshots.reserve(trace.partition_results.size());
  for (const auto &result : trace.partition_results) {
    mcpd3::DualDecompositionPartitionSnapshot snapshot;
    snapshot.partition_id = result.partition_id;
    snapshot.lower_bound = result.lower_bound;
    snapshot.regularization_budget = result.regularization_budget;
    snapshot.regularization_contribution =
        result.regularization_contribution;
    snapshot.regularization_anchor_sink_count =
        result.regularization_anchor_sink_count;
    snapshot.regularization_active_sink_count =
        result.regularization_active_sink_count;
    snapshot.local_labels.assign(result.full_labels.size(), 0);
    for (const auto &label : result.full_labels) {
      require(label.local_index >= 0 &&
                  static_cast<size_t>(label.local_index) <
                      snapshot.local_labels.size(),
              "trace full label local index out of range");
      snapshot.local_labels[static_cast<size_t>(label.local_index)] =
          label.label;
    }
    snapshots.push_back(std::move(snapshot));
  }
  std::sort(snapshots.begin(), snapshots.end(), [](const auto &lhs,
                                                   const auto &rhs) {
    return lhs.partition_id < rhs.partition_id;
  });
  return snapshots;
}

void requirePartitionSnapshotsEqual(
    const std::vector<mcpd3::DualDecompositionPartitionSnapshot> &actual,
    const std::vector<mcpd3::DualDecompositionPartitionSnapshot> &expected,
    const std::string &context) {
  require(actual.size() == expected.size(),
          context + ": partition snapshot count differs");
  for (size_t i = 0; i < actual.size(); ++i) {
    const auto &lhs = actual[i];
    const auto &rhs = expected[i];
    const std::string item_context =
        context + ": partition " + std::to_string(i);
    require(lhs.partition_id == rhs.partition_id,
            item_context + " id differs");
    require(lhs.lower_bound == rhs.lower_bound,
            item_context + " lower bound differs");
    require(lhs.regularization_budget == rhs.regularization_budget,
            item_context + " regularization budget differs");
    require(lhs.regularization_contribution ==
                rhs.regularization_contribution,
            item_context + " regularization contribution differs");
    require(lhs.regularization_anchor_sink_count ==
                rhs.regularization_anchor_sink_count,
            item_context + " regularization anchor count differs");
    require(lhs.regularization_active_sink_count ==
                rhs.regularization_active_sink_count,
            item_context + " regularization active count differs");
    require(lhs.local_labels == rhs.local_labels,
            item_context + " local labels differ");
  }
}

void partitionWorkerCoordinatorMatchesDualDecompositionRounds() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  auto package_source = makeTinyDualDecomposition();
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  for (size_t i = 0; i < package_source.getPartitionPackages().size(); ++i) {
    workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  }

  mcpd3::PartitionWorkerCoordinatorOptions coordinator_options;
  coordinator_options.min_step_size = 1;
  coordinator_options.use_momentum = false;
  mcpd3::PartitionWorkerCoordinator coordinator(
      package_source.getPartitionPackages(), std::move(workers),
      coordinator_options);
  auto reference = makeTinyDualDecomposition();
  requireConstraintSnapshotsEqual(
      coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
      "initial");

  mcpd3::Objective best_worker_lower_bound = 0;
  bool has_best_worker_lower_bound = false;
  mcpd3::PartitionWorkerRoundStats worker_stats;
  for (long round = 1; round <= 2; ++round) {
    const auto trace = coordinator.runRoundWithTrace(
        /*round_id=*/round, /*scale=*/100, /*step_size=*/100,
        /*regularization_strength=*/0);
    worker_stats = trace.stats;
    if (!has_best_worker_lower_bound ||
        worker_stats.lower_bound > best_worker_lower_bound) {
      best_worker_lower_bound = worker_stats.lower_bound;
      has_best_worker_lower_bound = true;
    }
    reference.runOptimizationScale(
        /*nstep=*/1, /*step_size=*/100, /*max_cycle_count=*/2,
        /*use_momentum=*/false);
    requireConstraintSnapshotsEqual(
        coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
        "round " + std::to_string(round));
    requirePartitionSnapshotsEqual(
        partitionSnapshotsFromTrace(trace), reference.getPartitionSnapshots(),
        "round " + std::to_string(round));
  }

  require(best_worker_lower_bound == reference.getBestLowerBoundRaw(),
          "coordinator best lower bound differs from DualDecomposition");
  require(worker_stats.disagreement_count ==
              reference.getLastDisagreementCount(),
          "coordinator disagreement count differs from DualDecomposition");
  require(worker_stats.disagreement_norm_sq ==
              reference.getLastDisagreementNormSq(),
          "coordinator disagreement norm differs from DualDecomposition");
}

mcpd3::MinCutGraph makeParityFixtureGraph() {
  mcpd3::MinCutGraph graph;
  graph.nnode = 8;
  graph.narc = 12;
  graph.arcs = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
                6, 7, 0, 4, 1, 5, 2, 6, 3, 7, 0, 7};
  graph.arc_capacities = {7,  3,  5,  11, 13, 2,  17, 19,
                          23, 5,  29, 7,  31, 37, 41, 3,
                          2,  43, 47, 11, 53, 13, 59, 17};
  graph.terminal_capacities = {19, -7, 11, -13, 17, -5, 23, -29};
  return graph;
}

mcpd3::DualDecompositionOptions makeParityDualOptions(bool use_momentum) {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 1000000;
  options.use_momentum = use_momentum;
  options.enable_group_stopping = false;
  options.exhaust_regularized_scale_iterations = true;
  return options;
}

mcpd3::PartitionWorkerCoordinatorOptions
makeParityCoordinatorOptions(const mcpd3::DualDecompositionOptions &dual) {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.objective_scale = dual.objective_scale;
  options.use_momentum = dual.use_momentum;
  options.enable_group_stopping = dual.enable_group_stopping;
  options.exhaust_regularized_scale_iterations =
      dual.exhaust_regularized_scale_iterations;
  options.regularization_budget_limit = dual.regularization_budget_limit;
  options.min_step_size = dual.min_step_size;
  options.regularization_scheme =
      dual.regularization_scheme ==
              mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON
          ? mcpd3::PartitionWorkerRegularizationScheme::SCALED_EPSILON
          : mcpd3::PartitionWorkerRegularizationScheme::NONE;
  return options;
}

void requireRoundTraceMatchesReference(
    const mcpd3::PartitionWorkerRoundTrace &trace,
    const mcpd3::DualDecomposition &reference, const std::string &context) {
  require(trace.stats.original_objective ==
              reference.getLastOriginalObjectiveRaw(),
          context + ": original objective differs");
  require(trace.stats.certified_lower_bound ==
              reference.getLastCertifiedLowerBoundRaw(),
          context + ": certified lower bound differs");
  require(trace.stats.lower_bound == reference.getLastCertifiedLowerBoundRaw(),
          context + ": lower bound differs");
  require(trace.stats.regularized_objective ==
              reference.getLastRegularizedObjectiveRaw(),
          context + ": regularized objective differs");
  require(trace.stats.regularization_budget ==
              reference.getLastRegularizationBudget(),
          context + ": regularization budget differs");
  require(trace.stats.regularization_contribution ==
              reference.getLastRegularizationContribution(),
          context + ": regularization contribution differs");
  require(trace.stats.regularization_anchor_sink_count ==
              reference.getLastRegularizationAnchorSinkCount(),
          context + ": regularization anchor count differs");
  require(trace.stats.regularization_active_sink_count ==
              reference.getLastRegularizationActiveSinkCount(),
          context + ": regularization active count differs");
  require(trace.stats.disagreement_count ==
              reference.getLastDisagreementCount(),
          context + ": disagreement count differs");
  require(trace.stats.disagreement_norm_sq ==
              reference.getLastDisagreementNormSq(),
          context + ": disagreement norm differs");
}

void partitionWorkerCoordinatorMatchesDualDecompositionParityFixture(
    bool use_momentum) {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  const auto options = makeParityDualOptions(use_momentum);
  auto package_source = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), options);

  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());

  mcpd3::PartitionWorkerCoordinator coordinator(
      package_source.getPartitionPackages(), std::move(workers),
      makeParityCoordinatorOptions(options));
  auto reference = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), options);

  requireConstraintSnapshotsEqual(
      coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
      "parity fixture initial");

  const std::vector<long> step_sizes = {1000, 100, 10, 10, 1, 1, 100};
  for (size_t round_index = 0; round_index < step_sizes.size();
       ++round_index) {
    const long step_size = step_sizes[round_index];
    const int regularization_strength =
        reference.regularizationStrengthForStepSize(step_size);
    const auto trace = coordinator.runRoundWithTrace(
        static_cast<long>(round_index + 1), options.objective_scale, step_size,
        regularization_strength);
    reference.runOptimizationScale(
        /*nstep=*/1, step_size, /*max_cycle_count=*/2, use_momentum);

    const std::string context =
        std::string(use_momentum ? "momentum" : "no momentum") +
        " parity fixture round " + std::to_string(round_index + 1) +
        " step " + std::to_string(step_size);
    requireRoundTraceMatchesReference(trace, reference, context);
    requireConstraintSnapshotsEqual(
        coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
        context);
    requirePartitionSnapshotsEqual(partitionSnapshotsFromTrace(trace),
                                   reference.getPartitionSnapshots(), context);
  }
}

void partitionWorkerCoordinatorMatchesDualDecompositionRegularizedRounds() {
  partitionWorkerCoordinatorMatchesDualDecompositionParityFixture(
      /*use_momentum=*/false);
  partitionWorkerCoordinatorMatchesDualDecompositionParityFixture(
      /*use_momentum=*/true);
}

void requireLocalWarmStatesEqual(
    const std::vector<mcpd3::PrimalDualMinCutSolver::WarmState> &actual,
    const std::vector<mcpd3::PrimalDualMinCutSolver::WarmState> &expected,
    const std::string &context);
constexpr bool fileBackedNativeDdSupported();

void partitionWorkerCoordinatorCapacityRefreshMatchesNativeState() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  auto options = makeParityDualOptions(/*use_momentum=*/true);
  options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  auto initial_graph = makeParityFixtureGraph();
  mcpd3::DualDecomposition reference(
      /*npartition=*/4, initial_graph, options);

  auto worker = std::make_unique<mcpd3::InProcessPartitionWorker>();
  auto *worker_ptr = worker.get();
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::move(worker));
  mcpd3::PartitionWorkerCoordinator coordinator(
      reference.getPartitionPackages(), std::move(workers),
      makeParityCoordinatorOptions(options));

  for (long round = 1; round <= 3; ++round) {
    const auto trace = coordinator.runRoundWithTrace(
        round, options.objective_scale, /*step_size=*/100,
        reference.regularizationStrengthForStepSize(100));
    reference.runOptimizationScale(/*nstep=*/1, /*step_size=*/100,
                                   /*max_cycle_count=*/2,
                                   /*use_momentum=*/true);
    requireRoundTraceMatchesReference(
        trace, reference, "capacity refresh prefix round " +
                              std::to_string(round));
  }

  auto replacement_graph = initial_graph;
  for (auto &capacity : replacement_graph.arc_capacities) {
    capacity *= 2;
  }
  for (auto &capacity : replacement_graph.terminal_capacities) {
    capacity *= 2;
  }

  auto package_options = options;
  package_options.construct_solvers = false;
  mcpd3::DualDecomposition replacement_packages(
      /*npartition=*/4, replacement_graph, package_options);
  std::vector<mcpd3::PartitionCapacityUpdate> updates;
  for (const auto &package : replacement_packages.getPartitionPackages()) {
    updates.push_back(mcpd3::PartitionCapacityUpdate{
        package.partition_id,
        package.arc_capacities.copyToVector(),
        package.terminal_capacities.copyToVector(),
        /*preserve_flow_state=*/true,
        /*flow_scale_numerator=*/2,
        /*flow_scale_denominator=*/1});
  }

  reference.replaceProblemCapacities(
      replacement_graph.arc_capacities,
      replacement_graph.terminal_capacities,
      /*preserve_alpha_state=*/true, /*preserve_flow_state=*/true,
      /*flow_scale_numerator=*/2, /*flow_scale_denominator=*/1);
  coordinator.replacePartitionCapacities(
      updates, /*preserve_alpha_state=*/true);

  std::vector<mcpd3::PrimalDualMinCutSolver::WarmState> worker_states;
  for (const auto &package : replacement_packages.getPartitionPackages()) {
    worker_states.push_back(worker_ptr->warmState(package.partition_id));
  }
  requireLocalWarmStatesEqual(worker_states,
                              reference.captureLocalSolverWarmStates(),
                              "capacity refresh immediate state");
  requireConstraintSnapshotsEqual(
      coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
      "capacity refresh preserved alpha");

  for (long round = 4; round <= 8; ++round) {
    const long step_size = round < 7 ? 10 : 1;
    const int regularization_strength =
        reference.regularizationStrengthForStepSize(step_size);
    const auto trace = coordinator.runRoundWithTrace(
        round, options.objective_scale, step_size, regularization_strength);
    reference.runOptimizationScale(/*nstep=*/1, step_size,
                                   /*max_cycle_count=*/2,
                                   /*use_momentum=*/true);
    const std::string context =
        "capacity refresh suffix round " + std::to_string(round);
    requireRoundTraceMatchesReference(trace, reference, context);
    requireConstraintSnapshotsEqual(
        coordinator.getConstraintSnapshots(), reference.getConstraintSnapshots(),
        context);
    requirePartitionSnapshotsEqual(partitionSnapshotsFromTrace(trace),
                                   reference.getPartitionSnapshots(), context);
  }

  auto same_capacity_updates = updates;
  for (auto &update : same_capacity_updates) {
    update.flow_scale_numerator = 1;
    update.flow_scale_denominator = 1;
  }
  coordinator.replacePartitionCapacities(
      same_capacity_updates, /*preserve_alpha_state=*/false);
  const auto reset_constraints = coordinator.getConstraintSnapshots();
  for (const auto &constraint : reset_constraints) {
    require(constraint.alpha == 0 && constraint.last_alpha == 0 &&
                constraint.alpha_momentum == 0,
            "capacity refresh alpha reset must clear all multiplier state");
  }

  requireThrows(
      [&] {
        auto missing = updates;
        missing.pop_back();
        coordinator.replacePartitionCapacities(missing);
      },
      "capacity refresh must reject a missing partition");
  requireThrows(
      [&] {
        auto duplicate = updates;
        duplicate.back().partition_id = duplicate.front().partition_id;
        coordinator.replacePartitionCapacities(duplicate);
      },
      "capacity refresh must reject a duplicate partition");
}

void coordinatorReconfiguresOnlyThePerSolveSchedule() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);
  auto dual_options = makeParityDualOptions(/*use_momentum=*/false);
  dual_options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  mcpd3::DualDecomposition reference(
      /*npartition=*/4, makeParityFixtureGraph(), dual_options);
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  auto options = makeParityCoordinatorOptions(dual_options);
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.initial_step_size = 1000;
  options.enable_group_stopping = false;
  options.promote_objective_scale_on_overbudget = false;
  options.regularization_scheme =
      mcpd3::PartitionWorkerRegularizationScheme::NONE;
  mcpd3::PartitionWorkerCoordinator coordinator(
      reference.getPartitionPackages(), std::move(workers), options);

  coordinator.configureOptimizationSchedule(
      /*num_optimization_scales=*/3, /*initial_step_size=*/100,
      /*exhaust_regularized_scale_iterations=*/true);
  const auto result = coordinator.solve();
  require(!result.scale_results.empty(),
          "reconfigured coordinator should execute a scale");
  require(result.scale_results.front().step_size == 100,
          "reconfigured coordinator should use the new initial step");
  require(coordinator.objectiveScale() == options.objective_scale,
          "schedule reconfiguration must not alter objective scale");
  requireThrows(
      [&] { coordinator.configureOptimizationSchedule(0, 1, false); },
      "schedule reconfiguration must reject zero scale count");
  requireThrows(
      [&] { coordinator.configureOptimizationSchedule(1, 0, false); },
      "schedule reconfiguration must reject zero initial step");
}

void inProcessWorkerFileBacksEveryPersistentSolverArray() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-worker-file-backed-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);

  mcpd3::PartitionPackage package;
  package.partition_id = 3;
  package.local_node_count = 2;
  package.arcs = {0, 1};
  package.arc_capacities = {7, 7};
  package.terminal_capacities = {11, -11};
  package.local_to_global = {4, 5};

  mcpd3::SolverStorageOptions storage;
  storage.mode = mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  storage.directory = scratch.string();
  if constexpr (!(mcpd3::solver_storage_mmap_compatible_v<mcpd3::Capacity> &&
                  mcpd3::solver_storage_mmap_compatible_v<mcpd3::NodeFlow> &&
                  mcpd3::solver_storage_mmap_compatible_v<
                      mcpd3::TerminalResidual> &&
                  mcpd3::solver_storage_mmap_compatible_v<mcpd3::Objective>)) {
    requireThrows(
        [&] {
          mcpd3::InProcessPartitionWorker unsupported(storage);
          unsupported.loadPartition(package);
        },
        "file-backed worker must reject nontrivial solver types");
    std::filesystem::remove_all(scratch);
    return;
  }

  mcpd3::InProcessPartitionWorker resident;
  mcpd3::InProcessPartitionWorker file_backed(storage);
  resident.loadPartition(package);
  file_backed.loadPartition(package);
  const auto diagnostics = file_backed.storageDiagnostics(3);
  require(diagnostics.topology_file_backed,
          "worker topology must be file-backed");
  require(diagnostics.capacity_file_backed,
          "worker capacities must be file-backed");
  require(diagnostics.flow_file_backed,
          "worker flow must be file-backed");
  require(diagnostics.labels_file_backed,
          "worker labels must be file-backed");
  require(diagnostics.bk_nodes_file_backed &&
              diagnostics.bk_arcs_file_backed,
          "worker BK state must be file-backed");
  require(diagnostics.bk_changed_arc_marks_file_backed,
          "worker BK changed-arc marks must be file-backed");
  require(diagnostics.file_backed_bytes > 0,
          "worker must report file-backed bytes");

  mcpd3::PartitionSolveRequest request;
  request.round_id = 1;
  request.partition_id = 3;
  request.return_full_labels = true;
  const auto expected = resident.solveRound(request);
  const auto actual = file_backed.solveRound(request);
  require(actual.lower_bound == expected.lower_bound,
          "file-backed worker lower bound must match resident");
  require(actual.full_labels.size() == expected.full_labels.size(),
          "file-backed worker label count must match resident");
  for (size_t index = 0; index < actual.full_labels.size(); ++index) {
    require(actual.full_labels[index].label == expected.full_labels[index].label,
            "file-backed worker labels must match resident");
  }
  std::filesystem::remove_all(scratch);
}

void requireFlowWarmStartsEqual(
    const mcpd3::DualDecomposition::FlowWarmStart &actual,
    const mcpd3::DualDecomposition::FlowWarmStart &expected,
    const std::string &context) {
  require(actual.constraints.size() == expected.constraints.size(),
          context + ": constraint count differs");
  requireConstraintSnapshotsEqual(actual.constraints, expected.constraints,
                                  context);
  require(actual.partitions.size() == expected.partitions.size(),
          context + ": partition count differs");
  for (size_t partition = 0; partition < actual.partitions.size(); ++partition) {
    const auto &actual_partition = actual.partitions[partition];
    const auto &expected_partition = expected.partitions[partition];
    const std::string partition_context =
        context + " partition " + std::to_string(partition);
    require(actual_partition.arcs == expected_partition.arcs,
            partition_context + ": topology differs");
    require(actual_partition.arc_capacities ==
                expected_partition.arc_capacities,
            partition_context + ": arc capacities differ");
    require(actual_partition.terminal_capacities ==
                expected_partition.terminal_capacities,
            partition_context + ": terminal capacities differ");
    require(actual_partition.v_flow == expected_partition.v_flow,
            partition_context + ": arc flow differs");
    require(actual_partition.d_flow == expected_partition.d_flow,
            partition_context + ": node flow differs");
    require(actual_partition.x == expected_partition.x,
            partition_context + ": labels differ");
  }
}

void requireLocalWarmStatesEqual(
    const std::vector<mcpd3::PrimalDualMinCutSolver::WarmState> &actual,
    const std::vector<mcpd3::PrimalDualMinCutSolver::WarmState> &expected,
    const std::string &context) {
  require(actual.size() == expected.size(),
          context + ": local warm-state count differs");
  for (size_t partition = 0; partition < actual.size(); ++partition) {
    const auto &a = actual[partition];
    const auto &e = expected[partition];
    const std::string p = context + " partition " + std::to_string(partition);
    require(a.v_flow == e.v_flow, p + ": arc flow differs");
    require(a.d_flow == e.d_flow, p + ": node flow differs");
    require(a.x == e.x, p + ": labels differ");
    require(a.is_first_iteration == e.is_first_iteration,
            p + ": first-iteration state differs");
    require(a.is_first_iteration_of_new_scale ==
                e.is_first_iteration_of_new_scale,
            p + ": new-scale state differs");
    require(a.has_solution == e.has_solution, p + ": solved state differs");
    require(a.mincut_value == e.mincut_value, p + ": mincut value differs");
    require(a.cached_lagrange_multipliers == e.cached_lagrange_multipliers,
            p + ": cached alpha differs");
    require(a.cached_last_lagrange_multipliers ==
                e.cached_last_lagrange_multipliers,
            p + ": cached last alpha differs");
    require(a.regularization_str == e.regularization_str,
            p + ": regularization strength differs");
    require(a.last_regularization_budget == e.last_regularization_budget,
            p + ": regularization budget differs");
    require(a.last_regularization_contribution ==
                e.last_regularization_contribution,
            p + ": regularization contribution differs");
    require(a.last_regularization_anchor_sink_count ==
                e.last_regularization_anchor_sink_count,
            p + ": regularization anchor count differs");
    require(a.last_regularization_active_sink_count ==
                e.last_regularization_active_sink_count,
            p + ": regularization active count differs");
    require(a.regularization_weights == e.regularization_weights,
            p + ": regularization weights differ");
    const auto &ag = a.maxflow_graph_state;
    const auto &eg = e.maxflow_graph_state;
    require(ag.node_num == eg.node_num && ag.arc_num == eg.arc_num,
            p + ": BK graph shape differs");
    require(ag.flow == eg.flow, p + ": BK flow differs");
    require(ag.maxflow_iteration == eg.maxflow_iteration,
            p + ": BK iteration differs");
    require(ag.time == eg.time, p + ": BK timestamp differs");
    require(ag.node_tr_caps == eg.node_tr_caps,
            p + ": BK terminal residuals differ");
    require(ag.node_parent_arc_indices == eg.node_parent_arc_indices,
            p + ": BK search-tree parents differ");
    require(ag.node_timestamps == eg.node_timestamps,
            p + ": BK node timestamps differ");
    require(ag.node_distances == eg.node_distances,
            p + ": BK tree distances differ");
    require(ag.node_is_sink == eg.node_is_sink,
            p + ": BK tree labels differ");
    require(ag.arc_residual_capacities == eg.arc_residual_capacities,
            p + ": BK residual capacities differ");
  }
}

constexpr bool fileBackedNativeDdSupported() {
  return mcpd3::solver_storage_mmap_compatible_v<mcpd3::Capacity> &&
         mcpd3::solver_storage_mmap_compatible_v<mcpd3::NodeFlow> &&
         mcpd3::solver_storage_mmap_compatible_v<mcpd3::TerminalResidual> &&
         mcpd3::solver_storage_mmap_compatible_v<mcpd3::Objective>;
}

void fileBackedNativeDdMatchesResidentStateExactly() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-native-file-backed-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);

  auto resident_options = makeParityDualOptions(/*use_momentum=*/true);
  resident_options.emit_partition_packages = false;
  resident_options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};

  auto file_options = resident_options;
  file_options.solver_storage.mode =
      mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  file_options.solver_storage.directory = scratch.string();

  if constexpr (!fileBackedNativeDdSupported()) {
    requireThrows(
        [&] {
          mcpd3::DualDecomposition unsupported(
              /*npartition=*/4, makeParityFixtureGraph(), file_options);
        },
        "file-backed DD must reject nontrivial capacity types");
    std::filesystem::remove_all(scratch);
    return;
  }

  auto resident = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), resident_options);
  auto file_backed = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), file_options);

  const auto diagnostics = file_backed.getLocalSolverStorageDiagnostics();
  require(diagnostics.size() == 4,
          "file-backed DD must report one storage record per partition");
  for (const auto &diagnostic : diagnostics) {
    require(diagnostic.topology_file_backed,
            "file-backed DD topology must use a file mapping");
    require(diagnostic.capacity_file_backed,
            "file-backed DD capacities must use file mappings");
    require(diagnostic.flow_file_backed,
            "file-backed DD primal-dual flow must use file mappings");
    require(diagnostic.labels_file_backed,
            "file-backed DD labels must use a file mapping");
    require(diagnostic.bk_nodes_file_backed,
            "file-backed DD BK nodes must use a file mapping");
    require(diagnostic.bk_arcs_file_backed,
            "file-backed DD BK residual arcs must use a file mapping");
    require(diagnostic.bk_changed_arc_marks_file_backed,
            "file-backed DD BK changed-arc marks must use a file mapping");
    require(diagnostic.file_backed_bytes > 0,
            "file-backed DD must report mapped state bytes");
  }

  const std::vector<long> step_sizes = {1000, 100, 10, 10, 1, 1, 100};
  for (size_t round = 0; round < step_sizes.size(); ++round) {
    resident.runOptimizationScale(/*nstep=*/1, step_sizes[round],
                                  /*max_cycle_count=*/2,
                                  /*use_momentum=*/true);
    file_backed.runOptimizationScale(/*nstep=*/1, step_sizes[round],
                                     /*max_cycle_count=*/2,
                                     /*use_momentum=*/true);
    const std::string context =
        "file-backed parity round " + std::to_string(round + 1);
    require(resident.getLastOriginalObjectiveRaw() ==
                file_backed.getLastOriginalObjectiveRaw(),
            context + ": original objective differs");
    require(resident.getLastCertifiedLowerBoundRaw() ==
                file_backed.getLastCertifiedLowerBoundRaw(),
            context + ": certified lower bound differs");
    require(resident.getLastRegularizedObjectiveRaw() ==
                file_backed.getLastRegularizedObjectiveRaw(),
            context + ": regularized objective differs");
    require(resident.getLastDisagreementCount() ==
                file_backed.getLastDisagreementCount(),
            context + ": disagreement count differs");
    require(resident.getLastDisagreementNormSq() ==
                file_backed.getLastDisagreementNormSq(),
            context + ": disagreement norm differs");
    require(resident.getLastRegularizationBudget() ==
                file_backed.getLastRegularizationBudget(),
            context + ": regularization budget differs");
    require(resident.getLastRegularizationContribution() ==
                file_backed.getLastRegularizationContribution(),
            context + ": regularization contribution differs");
    requireConstraintSnapshotsEqual(file_backed.getConstraintSnapshots(),
                                    resident.getConstraintSnapshots(), context);
    requirePartitionSnapshotsEqual(file_backed.getPartitionSnapshots(),
                                   resident.getPartitionSnapshots(), context);
    requireLocalWarmStatesEqual(file_backed.captureLocalSolverWarmStates(),
                                resident.captureLocalSolverWarmStates(),
                                context);
  }

  std::filesystem::remove_all(scratch);
}

void fileBackedNativeDdPreservesPersistentFlowExactly() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-native-file-backed-persistent-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);

  auto resident_options = makeParityDualOptions(/*use_momentum=*/true);
  resident_options.emit_partition_packages = false;
  resident_options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  resident_options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  auto file_options = resident_options;
  file_options.solver_storage.mode =
      mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  file_options.solver_storage.directory = scratch.string();

  if constexpr (!fileBackedNativeDdSupported()) {
    requireThrows(
        [&] {
          mcpd3::DualDecomposition unsupported(
              /*npartition=*/4, makeParityFixtureGraph(), file_options);
        },
        "persistent file-backed DD must reject nontrivial capacity types");
    std::filesystem::remove_all(scratch);
    return;
  }

  auto resident = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), resident_options);
  auto file_backed = mcpd3::DualDecomposition(
      /*npartition=*/4, makeParityFixtureGraph(), file_options);
  for (const long step : {1000L, 100L, 10L}) {
    resident.runOptimizationScale(1, step, 2, true);
    file_backed.runOptimizationScale(1, step, 2, true);
  }
  requireFlowWarmStartsEqual(file_backed.captureFlowWarmStart(),
                             resident.captureFlowWarmStart(),
                             "before persistent capacity replacement");

  const auto graph = makeParityFixtureGraph();
  std::vector<mcpd3::Capacity> replacement_arcs = graph.arc_capacities;
  std::vector<mcpd3::Capacity> replacement_terminals =
      graph.terminal_capacities;
  for (auto &capacity : replacement_arcs) {
    capacity *= 2;
  }
  for (auto &capacity : replacement_terminals) {
    capacity *= 2;
  }
  resident.replaceProblemCapacities(replacement_arcs, replacement_terminals,
                                    /*preserve_alpha_state=*/true,
                                    /*preserve_flow_state=*/true,
                                    /*flow_scale_numerator=*/2,
                                    /*flow_scale_denominator=*/1);
  file_backed.replaceProblemCapacities(
      replacement_arcs, replacement_terminals,
      /*preserve_alpha_state=*/true, /*preserve_flow_state=*/true,
      /*flow_scale_numerator=*/2, /*flow_scale_denominator=*/1);

  const auto resident_refresh = resident.getGlobalStorageDiagnostics();
  require(resident_refresh.last_capacity_refresh_arc_buffers_file_backed_partition_count ==
              0,
          "resident refresh buffers must remain resident");
  require(resident_refresh.last_capacity_refresh_terminal_buffers_file_backed_partition_count ==
              0,
          "resident terminal refresh buffers must remain resident");
  const auto mapped_refresh = file_backed.getGlobalStorageDiagnostics();
  require(mapped_refresh.last_capacity_refresh_arc_buffers_file_backed_partition_count ==
              4,
          "every file-backed local arc refresh buffer must be mapped");
  require(mapped_refresh.last_capacity_refresh_terminal_buffers_file_backed_partition_count ==
              4,
          "every file-backed local terminal refresh buffer must be mapped");
  for (const auto &diagnostic :
       file_backed.getLocalSolverStorageDiagnostics()) {
    require(diagnostic.last_capacity_replacement_flow_scratch_file_backed,
            "scaled flow replacement scratch must be file-backed");
  }

  for (const long step : {100L, 10L, 1L}) {
    resident.runOptimizationScale(1, step, 2, true);
    file_backed.runOptimizationScale(1, step, 2, true);
  }
  requireFlowWarmStartsEqual(file_backed.captureFlowWarmStart(),
                             resident.captureFlowWarmStart(),
                             "after persistent capacity replacement");
  requireConstraintSnapshotsEqual(file_backed.getConstraintSnapshots(),
                                  resident.getConstraintSnapshots(),
                                  "persistent file-backed constraints");
  requirePartitionSnapshotsEqual(file_backed.getPartitionSnapshots(),
                                 resident.getPartitionSnapshots(),
                                 "persistent file-backed partitions");
  requireLocalWarmStatesEqual(file_backed.captureLocalSolverWarmStates(),
                              resident.captureLocalSolverWarmStates(),
                              "persistent file-backed local state");

  std::filesystem::remove_all(scratch);
}

void fileBackedNativeDdRejectsInvalidStorage() {
  auto options = makeParityDualOptions(/*use_momentum=*/true);
  options.solver_storage.mode = mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  requireThrows(
      [&] {
        mcpd3::DualDecomposition invalid(
            /*npartition=*/4, makeParityFixtureGraph(), options);
      },
      "file-backed DD must require an explicit storage directory");

  options.solver_storage.directory =
      "/definitely/missing/mcpd3-native-file-backed-test";
  requireThrows(
      [&] {
        mcpd3::DualDecomposition invalid(
            /*npartition=*/4, makeParityFixtureGraph(), options);
      },
      "file-backed DD must reject a missing storage directory");
}

void nativeDdRecoversGlobalLabelsWithoutPartitionPackages() {
  auto package_options = makeParityDualOptions(/*use_momentum=*/true);
  package_options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  package_options.emit_partition_packages = true;
  package_options.retain_local_to_global_mapping = true;
  auto direct_options = package_options;
  direct_options.emit_partition_packages = false;

  mcpd3::DualDecomposition packaged(
      /*npartition=*/4, makeParityFixtureGraph(), package_options);
  mcpd3::DualDecomposition direct(
      /*npartition=*/4, makeParityFixtureGraph(), direct_options);
  packaged.solve();
  direct.solve();

  require(direct.getLastDisagreementCount() == 0,
          "package-free native DD must solve to agreement");
  const auto packaged_labels = packaged.getAgreedGlobalLabels();
  require(direct.getAgreedGlobalLabels() == packaged_labels,
          "package-free label recovery must match package-backed recovery");
  std::vector<int> direct_labels(packaged_labels.size(), -1);
  direct.copyAgreedGlobalLabels(direct_labels.data(), direct_labels.size());
  require(direct_labels == packaged_labels,
          "bounded global label copy must match vector recovery");
  requireThrows(
      [&] {
        direct.copyAgreedGlobalLabels(direct_labels.data(),
                                     direct_labels.size() - 1);
      },
      "global label copy must reject a mismatched destination size");
  requireThrows(
      [&] {
        direct.copyAgreedGlobalLabels(nullptr, direct_labels.size());
      },
      "global label copy must reject a null nonempty destination");
  requireThrows([&] { (void)direct.getPartitionPackages(); },
                "disabled package export must remain inaccessible");
}

void fileBackedNativeDdBacksPersistentGlobalMetadata() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-native-global-file-backed-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);

  auto options = makeParityDualOptions(/*use_momentum=*/true);
  options.emit_partition_packages = false;
  options.track_primal_upper_bound = true;
#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
  options.track_arc_flow_updates = true;
#endif
  options.materialize_all_partition_nodes = true;
  options.retain_local_to_global_mapping = true;
  options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  options.solver_storage.mode =
      mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  options.solver_storage.directory = scratch.string();

  if constexpr (!fileBackedNativeDdSupported()) {
    requireThrows(
        [&] {
          mcpd3::DualDecomposition unsupported(
              /*npartition=*/4, makeParityFixtureGraph(), options);
        },
        "file-backed global metadata must reject nontrivial capacity types");
    std::filesystem::remove_all(scratch);
    return;
  }

  mcpd3::DualDecomposition decomposition(
      /*npartition=*/4, makeParityFixtureGraph(), options);
  const auto diagnostics = decomposition.getGlobalStorageDiagnostics();
  require(diagnostics.original_topology_file_backed,
          "retained global topology must use file-backed storage");
  require(diagnostics.original_capacities_file_backed,
          "retained global capacities must use file-backed storage");
  require(diagnostics.partition_labels_file_backed,
          "partition labels must use file-backed storage");
  require(diagnostics.partition_validation_was_file_backed,
          "partition validation node markers must use file-backed storage");
  require(diagnostics.arc_locations_file_backed,
          "arc locations must use file-backed storage");
  require(diagnostics.terminal_locations_file_backed,
          "terminal locations must use file-backed storage");
  require(diagnostics.local_to_global_file_backed_partition_count == 4,
          "every retained local-to-global map must be file-backed");
  require(diagnostics.global_to_local_file_backed_partition_count == 4,
          "every retained construction lookup must be file-backed");
  require(diagnostics.local_graph_construction_file_backed_partition_count ==
              4,
          "every local graph must be constructed directly in file-backed "
          "storage");
  require(diagnostics.file_backed_bytes > 0,
          "global metadata must report mapped bytes");

  decomposition.solve();
  mcpd3::SolverArray<int> recovered_labels(
      8, options.solver_storage, "direct_global_label_recovery");
  decomposition.copyAgreedGlobalLabels(recovered_labels.data(),
                                       recovered_labels.size());
  require(recovered_labels.isFileBacked() &&
              recovered_labels.equals(decomposition.getAgreedGlobalLabels()),
          "global labels must recover directly into file-backed storage");
  std::filesystem::remove_all(scratch);
}

void fileBackedBasicPartitioningGeneratesLabelsWithoutResidentStaging() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-native-partition-label-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);
  ::setenv("MCPD3_PARTITIONER", "basic", 1);

  auto generated_options = makeParityDualOptions(/*use_momentum=*/true);
  generated_options.emit_partition_packages = false;
  generated_options.materialize_all_partition_nodes = true;
  generated_options.solver_storage.mode =
      mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  generated_options.solver_storage.directory = scratch.string();

  if constexpr (!fileBackedNativeDdSupported()) {
    requireThrows(
        [&] {
          mcpd3::DualDecomposition unsupported(
              /*npartition=*/4, makeParityFixtureGraph(), generated_options);
        },
        "file-backed partition labels must reject nontrivial capacity types");
    std::filesystem::remove_all(scratch);
    return;
  }

  mcpd3::DualDecomposition generated(
      /*npartition=*/4, makeParityFixtureGraph(), generated_options);
  const auto generated_diagnostics = generated.getGlobalStorageDiagnostics();
  require(generated_diagnostics.partition_labels_generated_in_backing_store,
          "basic partition labels must be generated directly in mmap storage");
  require(generated_diagnostics.partition_label_generation_was_file_backed,
          "basic partition label generation must use a file-backed mapping");

  auto supplied_options = generated_options;
  supplied_options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  mcpd3::DualDecomposition supplied(
      /*npartition=*/4, makeParityFixtureGraph(), supplied_options);
  require(!supplied.getGlobalStorageDiagnostics()
               .partition_labels_generated_in_backing_store,
          "caller-supplied partition labels must remain distinguishable from "
          "direct generation");

  std::filesystem::remove_all(scratch);
}

void packageOnlyFileBackedDdTransfersMappedPayloadWithoutCopies() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-package-payload-file-backed-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);

  auto options = makeParityDualOptions(/*use_momentum=*/true);
  options.construct_solvers = false;
  options.emit_partition_packages = true;
  options.materialize_all_partition_nodes = true;
  options.partition_labels = {0, 0, 1, 1, 2, 2, 3, 3};
  options.solver_storage.mode = mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  options.solver_storage.directory = scratch.string();

  if constexpr (!fileBackedNativeDdSupported()) {
    requireThrows(
        [&] {
          mcpd3::DualDecomposition unsupported(
              /*npartition=*/4, makeParityFixtureGraph(), options);
        },
        "file-backed package payloads must reject nontrivial capacity types");
    std::filesystem::remove_all(scratch);
    return;
  }

  mcpd3::DualDecomposition decomposition(
      /*npartition=*/4, makeParityFixtureGraph(), options);
  size_t nonempty_package_count = 0;
  for (const auto &package : decomposition.getPartitionPackages()) {
    if (package.local_node_count == 0) {
      continue;
    }
    ++nonempty_package_count;
    require(package.arcs.isFileBacked(),
            "package topology must remain file-backed");
    require(package.arc_capacities.isFileBacked(),
            "package arc capacities must remain file-backed");
    require(package.terminal_capacities.isFileBacked(),
            "package terminal capacities must remain file-backed");
    require(package.local_to_global.isFileBacked(),
            "package local-to-global mapping must remain file-backed");
    require(package.fileBackedBytes() > 0,
            "package must report its mapped payload bytes");
  }
  require(nonempty_package_count == 4,
          "fixture must export four nonempty mapped packages");
  require(decomposition.getGlobalStorageDiagnostics()
              .partition_package_zero_copy_transfer_count == 4,
          "package-only export must move each local mapped graph into its "
          "package");

  std::filesystem::remove_all(scratch);
}

void directedStreamingDimacsMatchesGeneralReaderValue() {
  const std::string path = "/tmp/mcpd3-directed-streaming-dimacs-test.max";
  {
    std::ofstream out(path);
    out << "c directed streaming reader test\n";
    out << "p max 4 4\n";
    out << "n 1 s\n";
    out << "n 4 t\n";
    out << "a 1 2 5\n";
    out << "a 2 3 7\n";
    out << "a 3 2 11\n";
    out << "a 3 4 17\n";
  }

  auto general = mcpd3::read_dimacs(path);
  auto directed_streaming = mcpd3::read_dimacs_directed_streaming(path);
  require(general.nnode == directed_streaming.nnode,
          "directed streaming reader should preserve node count");
  require(general.terminal_capacities ==
              directed_streaming.terminal_capacities,
          "directed streaming reader should preserve terminal capacities");

  mcpd3::PrimalDualMinCutSolver general_solver(std::move(general));
  mcpd3::PrimalDualMinCutSolver directed_streaming_solver(
      std::move(directed_streaming));
  require(general_solver.maxflow() == directed_streaming_solver.maxflow(),
          "directed streaming reader should preserve maxflow value");
  std::remove(path.c_str());
}

void dualDecompositionRegularizationSchemeControlsLowScaleStrength() {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10000;

  mcpd3::DualDecomposition scaled_epsilon(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  require(scaled_epsilon.regularizationStrengthForStepSize(100) == 0,
          "high scales should not use scaled epsilon regularization");
  require(scaled_epsilon.regularizationStrengthForStepSize(10) == 10,
          "step 10 should use scaled epsilon regularization");
  require(scaled_epsilon.regularizationStrengthForStepSize(1) == 1,
          "step 1 should use scaled epsilon regularization");

  options.scaled_epsilon_max_step_size = 12;
  mcpd3::DualDecomposition earlier_scaled_epsilon(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  require(earlier_scaled_epsilon.regularizationStrengthForStepSize(13) == 0,
          "scaled epsilon should remain disabled above its configured step");
  require(earlier_scaled_epsilon.regularizationStrengthForStepSize(12) == 12,
          "scaled epsilon should activate at its configured step");

  options.scaled_epsilon_strength_cap = 3;
  mcpd3::DualDecomposition capped_scaled_epsilon(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  require(capped_scaled_epsilon.regularizationStrengthForStepSize(12) == 3,
          "scaled epsilon should honor its configured strength cap");
  require(capped_scaled_epsilon.regularizationStrengthForStepSize(2) == 2,
          "a strength cap should not increase a smaller DD step");

  options.scaled_epsilon_strength_cap = -1;
  requireThrows(
      [&] {
        mcpd3::DualDecomposition invalid(
            /*npartition=*/2,
            /*nnode=*/2,
            /*narc=*/1,
            /*arcs=*/std::vector<int>{0, 1},
            /*arc_capacities=*/std::vector<int>{3, 5},
            /*terminal_capacities=*/std::vector<int>{2, -4}, options);
      },
      "scaled epsilon strength cap should be nonnegative");

  options.scaled_epsilon_strength_cap = 0;
  options.scaled_epsilon_max_step_size = 0;
  requireThrows(
      [&] {
        mcpd3::DualDecomposition invalid(
            /*npartition=*/2,
            /*nnode=*/2,
            /*narc=*/1,
            /*arcs=*/std::vector<int>{0, 1},
            /*arc_capacities=*/std::vector<int>{3, 5},
            /*terminal_capacities=*/std::vector<int>{2, -4}, options);
      },
      "scaled epsilon maximum step should be positive");

  options.scaled_epsilon_max_step_size = 10;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  mcpd3::DualDecomposition none(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  require(none.regularizationStrengthForStepSize(10) == 0,
          "NONE scheme should disable low-scale local regularization");
}

void disagreementPlateauTrackerRequiresAFullFlatWindow() {
  requireThrows(
      [] { mcpd3::DisagreementPlateauRegularizationTracker invalid(0); },
      "plateau tracker should reject nonpositive patience");
  mcpd3::DisagreementPlateauRegularizationTracker tracker(/*patience=*/2);
  requireThrows(
      [&] { (void)tracker.observe(/*iteration=*/-1, /*disagreement_count=*/5); },
      "plateau tracker should reject a negative iteration");
  requireThrows(
      [&] { (void)tracker.observe(/*iteration=*/0, /*disagreement_count=*/-1); },
      "plateau tracker should reject a negative disagreement count");
  require(!tracker.observe(/*iteration=*/0, /*disagreement_count=*/5),
          "first disagreement observation should establish the baseline");
  require(!tracker.observe(/*iteration=*/1, /*disagreement_count=*/5),
          "one flat iteration should not exhaust patience two");
  require(tracker.observe(/*iteration=*/2, /*disagreement_count=*/5),
          "two flat iterations should activate plateau regularization");
  require(tracker.active(),
          "plateau tracker should remain active after activation");
  require(!tracker.observe(/*iteration=*/3, /*disagreement_count=*/4),
          "disagreement progress should reset the pulse window");
  require(!tracker.observe(/*iteration=*/4, /*disagreement_count=*/4),
          "one flat iteration should not emit another regularization pulse");
  require(tracker.observe(/*iteration=*/5, /*disagreement_count=*/4),
          "each additional pulse should require another full flat window");
  require(!tracker.observe(/*iteration=*/6, /*disagreement_count=*/4),
          "regularization must not pulse on consecutive iterations");
}

void disagreementPlateauTrackerResetsOnProgressAndScaleReset() {
  mcpd3::DisagreementPlateauRegularizationTracker tracker(/*patience=*/2);
  require(!tracker.observe(/*iteration=*/0, /*disagreement_count=*/5),
          "first scale should establish its disagreement baseline");
  require(!tracker.observe(/*iteration=*/1, /*disagreement_count=*/4),
          "lower disagreement should reset plateau patience");
  require(tracker.bestDisagreementCount() == 4,
          "the tracker should retain the best disagreement count");
  require(tracker.iterationsSinceImprovement(/*iteration=*/1) == 0,
          "new disagreement progress should reset its patience age");
  require(!tracker.observe(/*iteration=*/2, /*disagreement_count=*/4),
          "one flat iteration after progress should not activate");
  require(tracker.iterationsSinceImprovement(/*iteration=*/2) == 1,
          "a flat iteration should advance the disagreement patience age");
  require(tracker.observe(/*iteration=*/3, /*disagreement_count=*/4),
          "a full flat window after progress should activate");

  tracker.reset();
  require(!tracker.active(),
          "a new optimization scale should start with activation disabled");
  require(!tracker.observe(/*iteration=*/0, /*disagreement_count=*/4),
          "a reset scale should establish a new disagreement baseline");
  require(!tracker.observe(/*iteration=*/1, /*disagreement_count=*/3),
          "new-scale disagreement progress should reset patience again");
  require(!tracker.observe(/*iteration=*/2, /*disagreement_count=*/3),
          "new-scale tracker should require its complete flat window");
  require(tracker.observe(/*iteration=*/3, /*disagreement_count=*/3),
          "new-scale tracker should activate after its own flat window");
}

void disagreementPlateauOptionsValidateAndUseUnitStrength() {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10000;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::
          DISAGREEMENT_PLATEAU_EPSILON;
  options.disagreement_patience = 0;
  requireThrows(
      [&] {
        mcpd3::DualDecomposition invalid(
            /*npartition=*/2, /*nnode=*/2, /*narc=*/1,
            /*arcs=*/std::vector<int>{0, 1},
            /*arc_capacities=*/std::vector<int>{3, 5},
            /*terminal_capacities=*/std::vector<int>{2, -4}, options);
      },
      "plateau mode should reject nonpositive disagreement patience");

  options.disagreement_patience = 2;
  mcpd3::DualDecomposition plateau(
      /*npartition=*/2, /*nnode=*/2, /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  require(plateau.regularizationStrengthForStepSize(10000) == 0,
          "plateau mode should start unregularized at high scales");
  require(plateau.plateauRegularizationStrength() == 1,
          "plateau mode should use cumulative unit regularization at every "
          "scale");
}

void disagreementPlateauModeActivatesOnARealDdPlateau() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 100;
  options.initial_step_size = 100;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 5;
  options.patience = 2;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::
          DISAGREEMENT_PLATEAU_EPSILON;
  options.disagreement_patience = 1;

  mcpd3::DualDecomposition plateau(
      /*npartition=*/2, /*nnode=*/2, /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{0, 0},
      /*terminal_capacities=*/std::vector<int>{-100, -10}, options);
  plateau.solve();

  require(plateau.getDisagreementPlateauActivationCount() == 1,
          "a persistent real DD disagreement should activate plateau "
          "regularization exactly once in one scale");
  require(plateau.getLastRegularizationBudget() > 0,
          "the solve after plateau activation should consume cumulative unit "
          "regularization budget");
  require(plateau.getLastRegularizationBudget() < plateau.getScale(),
          "the plateau fixture must retain a strictly sub-quantum budget");
}

void dualDecompositionRandomizesExportedInitialAlphas() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.randomize_initial_alphas = true;
  options.initial_alpha_random_radius = 9;
  options.initial_alpha_random_seed = 9;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);

  mcpd3::Lagrange randomized_alpha = 0;
  int endpoint_count = 0;
  for (const auto &package : dual_decomp.getPartitionPackages()) {
    for (const auto &endpoint : package.constraint_endpoints) {
      if (endpoint_count == 0) {
        randomized_alpha = endpoint.alpha;
      }
      require(endpoint.alpha == randomized_alpha,
              "random initial alpha should match across endpoints");
      require(endpoint.last_alpha == randomized_alpha,
              "random initial alpha should update last_alpha with alpha");
      require(endpoint.alpha >= -options.initial_alpha_random_radius &&
                  endpoint.alpha <= options.initial_alpha_random_radius,
              "random initial alpha should stay inside configured radius");
      ++endpoint_count;
    }
  }
  require(endpoint_count == 2,
          "tiny decomposition should export one pairwise constraint");
  require(randomized_alpha != 0,
          "seeded random initial alpha should perturb the initial multiplier");

  mcpd3::DualDecompositionOptions invalid_options = options;
  invalid_options.initial_alpha_random_radius = -1;
  bool threw = false;
  try {
    mcpd3::DualDecomposition invalid(
        /*npartition=*/2,
        /*nnode=*/2,
        /*narc=*/1,
        /*arcs=*/std::vector<int>{0, 1},
        /*arc_capacities=*/std::vector<int>{3, 5},
        /*terminal_capacities=*/std::vector<int>{2, -4}, invalid_options);
  } catch (const std::runtime_error &e) {
    threw =
        std::string(e.what()).find("initial alpha random radius") !=
        std::string::npos;
  }
  require(threw, "negative initial alpha random radius should be rejected");
}

void dualDecompositionObjectiveScaleIsIndependentOfStepSize() {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.initial_step_size = 100;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{3, 5},
      /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  dual_decomp.solve();
  require(dual_decomp.getScale() == 1,
          "objective reporting scale should not follow the DD step size");

  options.objective_scale = 7;
  mcpd3::DualDecomposition scaled_dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{21, 35},
      /*terminal_capacities=*/std::vector<int>{14, -28}, options);
  scaled_dual_decomp.solve();
  require(scaled_dual_decomp.getScale() == 7,
          "explicit objective scale should be preserved");

  options.objective_scale = 0;
  bool threw = false;
  try {
    mcpd3::DualDecomposition invalid_scale_dual_decomp(
        /*npartition=*/2,
        /*nnode=*/2,
        /*narc=*/1,
        /*arcs=*/std::vector<int>{0, 1},
        /*arc_capacities=*/std::vector<int>{3, 5},
        /*terminal_capacities=*/std::vector<int>{2, -4}, options);
  } catch (const std::runtime_error &e) {
    threw = std::string(e.what()).find("objective scale") !=
            std::string::npos;
  }
  require(threw, "nonpositive dual decomposition objective scale should fail");
}

void dualDecompositionPolyakStepUsesPrimalDualGap() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.initial_step_size = 10000;
  options.max_iteration_count = 100;
  options.max_total_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.step_policy =
      mcpd3::DualDecompositionStepPolicy::PolyakUpperBoundGap;
  options.polyak_step_theta = 1.0;

  std::vector<mcpd3::DualDecompositionIterationRecord> records;
  options.iteration_callback =
      [&](const mcpd3::DualDecompositionIterationRecord &record) {
        records.push_back(record);
      };

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/4, makeParityFixtureGraph(), options);
  bool rejected_negative_upper_bound = false;
  try {
    dual_decomp.setKnownPrimalUpperBoundRaw(mcpd3::Objective{-1});
  } catch (const std::runtime_error &) {
    rejected_negative_upper_bound = true;
  }
  require(rejected_negative_upper_bound,
          "a known primal upper bound must be nonnegative");
  dual_decomp.setKnownPrimalUpperBoundRaw(mcpd3::Objective{1000});
  dual_decomp.solve();

  require(records.size() == 1,
          "Polyak probe should report exactly one DD iteration");
  const auto &record = records.front();
  require(record.disagreement_norm_sq > 0.0,
          "Polyak fixture must produce a nonzero subgradient");
  require(record.positive_disagreement_count +
                  record.negative_disagreement_count ==
              record.disagreement_count,
          "signed disagreement counts must cover every disagreement");
  require(record.disagreement_entered_count == 0 &&
              record.disagreement_exited_count == 0 &&
              record.disagreement_flipped_count == 0,
          "the first DD round must not report disagreement-set churn");
  require(dual_decomp.hasCurrentUpperBound(),
          "Polyak step requires a feasible primal upper bound");
  const auto gap = dual_decomp.getCurrentUpperBoundRaw() -
                   record.certified_lower_bound_raw;
  require(gap > 0, "Polyak fixture must have a positive primal-dual gap");
  const long expected = std::clamp(
      static_cast<long>(std::llround(
          mcpd3::integer_to_double(gap) / record.disagreement_norm_sq)),
      options.min_step_size, options.initial_step_size);
  require(record.effective_step_size == expected,
          "Polyak step must equal the rounded gap over subgradient norm");
  require(record.effective_step_size < record.step_size,
          "Polyak fixture should reduce the scheduled alpha step");
}

void exactBoundaryCycleDetectorRejectsHashOnlyMatches() {
  mcpd3::ExactBoundaryCycleDetector detector(
      /*max_period=*/4, /*minimum_repetitions=*/2);

  const auto observe = [&](std::vector<std::uint8_t> labels,
                           std::vector<std::int8_t> diffs) {
    return detector.observe(
        mcpd3::BoundaryCycleSample{/*state_hash=*/17,
                                   /*labels=*/std::move(labels),
                                   /*diffs=*/std::move(diffs)});
  };

  require(!observe({0, 1}, {1}).has_value(),
          "one state cannot establish a cycle");
  require(!observe({1, 0}, {-1}).has_value(),
          "two distinct states cannot establish a repeated cycle");
  require(!observe({0, 1}, {1}).has_value(),
          "an incomplete second period cannot establish a cycle");
  require(!observe({1, 1}, {-1}).has_value(),
          "matching hashes must not hide an exact-label mismatch");
}

void exactBoundaryCycleDetectorFindsRepeatedNonconstantSequence() {
  mcpd3::ExactBoundaryCycleDetector detector(
      /*max_period=*/4, /*minimum_repetitions=*/2);
  const mcpd3::BoundaryCycleSample positive{
      /*state_hash=*/11, /*labels=*/{0, 1, 0}, /*diffs=*/{1, 0}};
  const mcpd3::BoundaryCycleSample negative{
      /*state_hash=*/29, /*labels=*/{1, 0, 1}, /*diffs=*/{-1, 0}};

  require(!detector.observe(positive).has_value(),
          "first positive state cannot establish a cycle");
  require(!detector.observe(negative).has_value(),
          "first negative state cannot establish a cycle");
  require(!detector.observe(positive).has_value(),
          "partial second period cannot establish a cycle");
  const auto cycle = detector.observe(negative);
  require(cycle.has_value(), "two exact periods must establish a cycle");
  require(cycle->period == 2, "alternating states must have period two");
  require(cycle->samples.size() == 2,
          "detected cycle must retain one complete period");
  require(cycle->samples[0].labels == positive.labels &&
              cycle->samples[1].labels == negative.labels,
          "detected cycle must preserve state order");
  require(cycle->nextSample(0).labels == positive.labels &&
              cycle->nextSample(1).labels == negative.labels &&
              cycle->nextSample(2).labels == positive.labels,
          "cycle replay must wrap in the observed order");
}

void exactBoundaryCycleDetectorRejectsConstantPlateaus() {
  mcpd3::ExactBoundaryCycleDetector detector(
      /*max_period=*/4, /*minimum_repetitions=*/2);
  const mcpd3::BoundaryCycleSample unchanged{
      /*state_hash=*/5, /*labels=*/{0, 1}, /*diffs=*/{1}};
  require(!detector.observe(unchanged).has_value(),
          "first plateau state cannot establish a cycle");
  require(!detector.observe(unchanged).has_value(),
          "a constant plateau must not be treated as a cycle");
}

void exactBoundaryCycleDetectorHonorsRequiredRepetitions() {
  mcpd3::ExactBoundaryCycleDetector detector(
      /*max_period=*/2, /*minimum_repetitions=*/3);
  const mcpd3::BoundaryCycleSample a{
      /*state_hash=*/1, /*labels=*/{0}, /*diffs=*/{1}};
  const mcpd3::BoundaryCycleSample b{
      /*state_hash=*/2, /*labels=*/{1}, /*diffs=*/{-1}};
  require(!detector.observe(a).has_value(), "first sample cannot cycle");
  require(!detector.observe(b).has_value(), "first period cannot cycle");
  require(!detector.observe(a).has_value(), "second period is incomplete");
  require(!detector.observe(b).has_value(),
          "two periods must not satisfy a three-period requirement");
  require(!detector.observe(a).has_value(), "third period is incomplete");
  require(detector.observe(b).has_value(),
          "three complete periods must satisfy the requirement");
}

void exactBoundaryCycleDetectorRejectsInvalidInputs() {
  requireThrows(
      [] {
        mcpd3::ExactBoundaryCycleDetector detector(
            /*max_period=*/1, /*minimum_repetitions=*/2);
      },
      "cycle detection must reject a period bound below two");
  requireThrows(
      [] {
        mcpd3::ExactBoundaryCycleDetector detector(
            /*max_period=*/2, /*minimum_repetitions=*/1);
      },
      "cycle detection must reject fewer than two repetitions");

  mcpd3::ExactBoundaryCycleDetector detector(
      /*max_period=*/2, /*minimum_repetitions=*/2);
  requireThrows(
      [&] {
        (void)detector.observe(
            mcpd3::BoundaryCycleSample{/*state_hash=*/0,
                                       /*labels=*/{},
                                       /*diffs=*/{1}});
      },
      "cycle detection must reject missing exact labels");
  requireThrows(
      [&] {
        (void)detector.observe(
            mcpd3::BoundaryCycleSample{/*state_hash=*/0,
                                       /*labels=*/{1},
                                       /*diffs=*/{}});
      },
      "cycle detection must reject missing replay updates");

  const mcpd3::ExactBoundaryCycle invalid_cycle;
  requireThrows([&] { (void)invalid_cycle.nextSample(0); },
                "cycle replay must reject an empty cycle");
}

void speculativeCycleProbeAcceptanceCoversToleranceBoundary() {
  require(mcpd3::accept_speculative_cycle_probe(
              /*candidate_lower_bound=*/101,
              /*baseline_lower_bound=*/100,
              /*rollback_tolerance=*/0),
          "an improving speculative lower bound must be accepted");
  require(mcpd3::accept_speculative_cycle_probe(100, 100, 0),
          "an equal speculative lower bound must be accepted");
  require(!mcpd3::accept_speculative_cycle_probe(99, 100, 0),
          "a worse speculative lower bound must be rejected");
  require(mcpd3::accept_speculative_cycle_probe(98, 100, 2),
          "a candidate at the rollback tolerance must be accepted");
  require(!mcpd3::accept_speculative_cycle_probe(97, 100, 2),
          "a candidate outside the rollback tolerance must be rejected");
  requireThrows(
      [] { (void)mcpd3::accept_speculative_cycle_probe(1, 1, -1); },
      "negative speculative rollback tolerance must be rejected");
}

void speculativeCycleReplayRetainsTheLastSolvedAlpha() {
  mcpd3::Lagrange last_alpha = mcpd3::lagrange_from_integer(3);
  mcpd3::record_lagrange_update_origin(
      mcpd3::lagrange_from_integer(7), last_alpha,
      /*advance_last_alpha=*/true);
  require(last_alpha == mcpd3::lagrange_from_integer(7),
          "a real DD update must record its starting alpha");

  mcpd3::record_lagrange_update_origin(
      mcpd3::lagrange_from_integer(11), last_alpha,
      /*advance_last_alpha=*/false);
  require(last_alpha == mcpd3::lagrange_from_integer(7),
          "virtual replay must retain the alpha used by the last exact solve");
}

#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
void dualDecompositionReportsPerPartitionMaxflowWork() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 2;
  options.max_iteration_count = 10;
  options.max_total_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.initial_step_size = 10;
  options.patience = 9;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.track_maxflow_work_telemetry = true;

  std::vector<mcpd3::DualDecompositionIterationRecord> records;
  options.iteration_callback =
      [&](const mcpd3::DualDecompositionIterationRecord &record) {
        records.push_back(record);
      };

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/4, makeParityFixtureGraph(), options);
  dual_decomp.solve();

  require(records.size() == 1,
          "work telemetry probe should report one DD round");
  const auto &record = records.front();
  require(record.partition_solves.size() == 4,
          "work telemetry must report every local partition solve");
  mcpd3::Objective local_lower_bound_sum = 0;
  for (size_t index = 0; index < record.partition_solves.size(); ++index) {
    const auto &partition = record.partition_solves[index];
    require(partition.partition_index == index,
            "partition work telemetry must retain stable partition order");
    require(partition.work.enabled,
            "partition work telemetry must identify instrumented solves");
    require(!partition.work.maxflow.reused_trees,
            "the first local partition solve must be cold");
    require(partition.solve_microseconds >=
                partition.work.maxflow_microseconds,
            "partition wall must include its maxflow work");
    require(partition.work.tree_nodes_changed ==
                partition.work.maxflow.changed_tree_nodes,
            "solver and BK changed-tree counts must agree");
    local_lower_bound_sum = mcpd3::checked_add(
        local_lower_bound_sum, partition.local_lower_bound_raw,
        "test local lower-bound sum overflow");
  }
  require(local_lower_bound_sum == record.regularized_objective_raw,
          "per-partition objectives must sum to the DD round objective");
}
#endif

void coherentBoundaryDirectionRequiresDominantPartitionMode() {
  require(mcpd3::coherent_boundary_direction(
              /*positive_count=*/80, /*negative_count=*/0,
              /*constraint_count=*/100, /*dominance_fraction=*/0.8) == 1,
          "a dominant positive partition mode should be detected");
  require(mcpd3::coherent_boundary_direction(
              /*positive_count=*/0, /*negative_count=*/80,
              /*constraint_count=*/100, /*dominance_fraction=*/0.8) == -1,
          "a dominant negative partition mode should be detected");
  require(mcpd3::coherent_boundary_direction(
              /*positive_count=*/79, /*negative_count=*/0,
              /*constraint_count=*/100, /*dominance_fraction=*/0.8) == 0,
          "a minority partition mode must not be classified as coherent");
  require(mcpd3::coherent_boundary_direction(
              /*positive_count=*/0, /*negative_count=*/0,
              /*constraint_count=*/0, /*dominance_fraction=*/0.8) == 0,
          "an empty boundary must not have a coherent direction");
  require(mcpd3::coherent_boundary_bisection_step(
              /*current_step=*/100, /*minimum_step=*/1,
              /*previous_direction=*/-1, /*current_direction=*/1) == 50,
          "a coherent direction reversal must bisect the alpha step");
  require(mcpd3::coherent_boundary_bisection_step(
              /*current_step=*/100, /*minimum_step=*/1,
              /*previous_direction=*/1, /*current_direction=*/1) == 100,
          "a persistent coherent direction must preserve the alpha step");
  require(mcpd3::coherent_boundary_bisection_step(
              /*current_step=*/1, /*minimum_step=*/1,
              /*previous_direction=*/-1, /*current_direction=*/1) == 1,
          "coherent bisection must respect the minimum alpha step");

  for (const double invalid_fraction : {0.5, 1.01}) {
    mcpd3::DualDecompositionOptions options;
    options.coherent_boundary_dominance_fraction = invalid_fraction;
    bool rejected = false;
    try {
      mcpd3::DualDecomposition invalid(
          /*npartition=*/1, /*nnode=*/1, /*narc=*/0,
          /*arcs=*/std::vector<int>{},
          /*arc_capacities=*/std::vector<int>{},
          /*terminal_capacities=*/std::vector<int>{0}, options);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    require(rejected,
            "coherent boundary dominance must stay in (0.5, 1]");
  }
}

void disagreementTransitionClassificationCoversEveryStateChange() {
  require(mcpd3::classify_disagreement_transition(0, 1) ==
              mcpd3::DisagreementTransition::Entered,
          "zero-to-active disagreement should enter");
  require(mcpd3::classify_disagreement_transition(-1, 0) ==
              mcpd3::DisagreementTransition::Exited,
          "active-to-zero disagreement should exit");
  require(mcpd3::classify_disagreement_transition(-1, 1) ==
              mcpd3::DisagreementTransition::Flipped,
          "opposite active signs should flip");
  require(mcpd3::classify_disagreement_transition(1, 1) ==
              mcpd3::DisagreementTransition::Unchanged,
          "a persistent sign should remain unchanged");
  require(mcpd3::classify_disagreement_transition(0, 0) ==
              mcpd3::DisagreementTransition::Unchanged,
          "persistent agreement should remain unchanged");
}

void dualDecompositionProbeHonorsTotalIterationBudget() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10;
  options.initial_step_size = 100;
  options.num_optimization_scales = 3;
  options.max_iteration_count = 100;
  options.max_total_iteration_count = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.max_objective_scale_promotions = 4;

  std::vector<mcpd3::DualDecompositionIterationRecord> records;
  options.iteration_callback =
      [&](const mcpd3::DualDecompositionIterationRecord &record) {
        records.push_back(record);
      };

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{30, 50},
      /*terminal_capacities=*/std::vector<int>{20, -40}, options);

  dual_decomp.solve();

  require(dual_decomp.getTotalOptimizationIterations() == 1,
          "probe should stop at the global DD iteration budget");
  require(dual_decomp.getObjectiveScalePromotionCount() == 0,
          "an exhausted probe budget must not promote the objective");
  require(records.size() == 1,
          "probe callback should report every completed DD iteration");
  const auto &record = records.front();
  require(record.total_iteration == 1,
          "probe callback total iteration mismatch");
  require(record.scale_iteration == 1,
          "probe callback scale iteration mismatch");
  require(record.objective_scale == 10,
          "probe callback objective scale mismatch");
  require(record.step_size == 100,
          "probe callback scheduled step mismatch");
  require(record.best_certified_lower_bound_raw ==
              dual_decomp.getBestCertifiedLowerBoundRaw(),
          "probe callback should expose the best certified lower bound");
  require(record.disagreement_count ==
              dual_decomp.getLastDisagreementCount(),
          "probe callback should expose the disagreement count");
}

void dualDecompositionRejectsNegativeTotalIterationBudget() {
  mcpd3::DualDecompositionOptions options;
  options.max_total_iteration_count = -1;
  bool threw = false;
  try {
    mcpd3::DualDecomposition invalid(
        /*npartition=*/1,
        /*nnode=*/1,
        /*narc=*/0,
        /*arcs=*/std::vector<int>{},
        /*arc_capacities=*/std::vector<int>{},
        /*terminal_capacities=*/std::vector<int>{0}, options);
  } catch (const std::runtime_error &e) {
    threw = std::string(e.what()).find("total iteration") !=
            std::string::npos;
  }
  require(threw, "negative total DD iteration budget should be rejected");
}

void dualDecompositionUsesExplicitPartitionLabels() {
  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.partition_labels = {0, 1, 0, 1};

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/4,
      /*narc=*/4,
      /*arcs=*/std::vector<int>{0, 1, 1, 2, 2, 3, 3, 0},
      /*arc_capacities=*/std::vector<int>{1, 1, 1, 1, 1, 1, 1, 1},
      /*terminal_capacities=*/std::vector<int>{1, 0, 0, -1}, options);

  require(dual_decomp.getPartitionLabels() == options.partition_labels,
          "explicit partition labels should bypass graph partitioning");
}

void dualDecompositionRejectsInvalidExplicitPartitionLabels() {
  const auto rejects = [](std::vector<int> labels, int partition_count,
                          const std::string &message) {
    mcpd3::DualDecompositionOptions options;
    options.partition_labels = std::move(labels);
    bool threw = false;
    try {
      mcpd3::DualDecomposition invalid(
          partition_count,
          /*nnode=*/3,
          /*narc=*/2,
          /*arcs=*/std::vector<int>{0, 1, 1, 2},
          /*arc_capacities=*/std::vector<int>{1, 1, 1, 1},
          /*terminal_capacities=*/std::vector<int>{1, 0, -1}, options);
    } catch (const std::runtime_error &error) {
      threw = std::string(error.what()).find("partition label") !=
              std::string::npos;
    }
    require(threw, message);
  };

  rejects({0, 1}, 2,
          "explicit partition labels must cover every graph node");
  rejects({0, 2, 1}, 2,
          "explicit partition labels must be in the partition range");
  rejects({0, 0, 0}, 2,
          "explicit partition labels must populate every partition");

  mcpd3::DualDecompositionOptions options;
  options.partition_labels = {0, 1, 1};
  options.partition_edge_weights = {1, 1};
  bool threw = false;
  try {
    mcpd3::DualDecomposition invalid(
        /*npartition=*/2,
        /*nnode=*/3,
        /*narc=*/2,
        /*arcs=*/std::vector<int>{0, 1, 1, 2},
        /*arc_capacities=*/std::vector<int>{1, 1, 1, 1},
        /*terminal_capacities=*/std::vector<int>{1, 0, -1}, options);
  } catch (const std::runtime_error &error) {
    threw = std::string(error.what()).find("mutually exclusive") !=
            std::string::npos;
  }
  require(threw,
          "explicit labels must not silently ignore partition edge weights");
}

void dualDecompositionPromotesObjectiveScaleOnOverBudget() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10;
  options.initial_step_size = 10;
  options.num_optimization_scales = 3;
  options.max_iteration_count = 20;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.max_objective_scale_promotions = 2;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{0, 0},
      /*terminal_capacities=*/std::vector<int>{-100, -10}, options);

  dual_decomp.solve();

  require(dual_decomp.getObjectiveScalePromotionCount() == 1,
          "over-budget scaled epsilon solve should promote objective scale");
  require(dual_decomp.getScale() == 100,
          "promotion should increase objective scale by one decade");
  require(dual_decomp.getBestLowerBoundRaw() == 0,
          "over-budget regularized lower bound should not be accepted: got " +
              mcpd3::integer_to_string(
                  dual_decomp.getBestLowerBoundRaw()));
  require(dual_decomp.getLastRegularizationBudget() <
              dual_decomp.getScale(),
          "promoted solve should finish under the active budget limit");
  require(dual_decomp.getLastDisagreementCount() == 0,
          "promoted solve should preserve progress and reach agreement");
}

void dualDecompositionPromotesAfterUnitScaleExhaustion() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.max_objective_scale_promotions = 4;

  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{30, 50},
      /*terminal_capacities=*/std::vector<int>{20, -40}, options);

  dual_decomp.solve();

  require(dual_decomp.getObjectiveScalePromotionCount() > 0,
          "native unit-scale exhaustion should promote the objective");
  require(dual_decomp.getLastDisagreementCount() == 0,
          "native promoted schedule should eventually reach agreement");
  require(dual_decomp.getConfiguredInitialStepSize() == dual_decomp.getScale(),
          "native promotion should retain its effective restart step");
}

void dualDecompositionCanRetryUnitStepWithoutMomentum() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  mcpd3::DualDecompositionOptions options;
  options.track_primal_upper_bound = false;
  options.verbose = false;
  options.thread_count = 1;
  options.objective_scale = 10;
  options.initial_step_size = 1;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = true;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.max_objective_scale_promotions = 0;

  mcpd3::DualDecomposition baseline(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{30, 50},
      /*terminal_capacities=*/std::vector<int>{20, -40}, options);

  baseline.solve();

  require(baseline.getUnitStepNoMomentumRetryCount() == 0,
          "unit-step cleanup should remain disabled by default");
  require(baseline.getTotalOptimizationIterations() == 1,
          "disabled unit-step cleanup should retain the original schedule");

  options.retry_unit_step_without_momentum = true;
  mcpd3::DualDecomposition dual_decomp(
      /*npartition=*/2,
      /*nnode=*/2,
      /*narc=*/1,
      /*arcs=*/std::vector<int>{0, 1},
      /*arc_capacities=*/std::vector<int>{30, 50},
      /*terminal_capacities=*/std::vector<int>{20, -40}, options);

  dual_decomp.solve();

  require(dual_decomp.getUnitStepNoMomentumRetryCount() == 1,
          "unit-step cleanup should run once after momentum fails at scale 1");
  require(dual_decomp.getTotalOptimizationIterations() == 2,
          "unit-step cleanup should add exactly one scale attempt");
}

struct ScriptedRound {
  long lower_bound = 0;
  int label = 0;
  long regularization_budget = 0;
  long regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
};

class ScriptedPartitionWorker final : public mcpd3::PartitionWorker {
public:
  explicit ScriptedPartitionWorker(std::deque<ScriptedRound> script)
      : script_(std::move(script)) {}

  void loadPartition(const mcpd3::PartitionPackage &package) override {
    packages_[package.partition_id] = package;
  }

  mcpd3::PartitionSolveResult solveRound(
      const mcpd3::PartitionSolveRequest &request) override {
    require(!script_.empty(), "scripted worker was called too many times");
    const auto &package = packageForRequest(request);
    requests_.push_back(request);
    const auto round = script_.front();
    script_.pop_front();

    mcpd3::PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = package.partition_id;
    result.lower_bound = round.lower_bound;
    result.regularization_budget = round.regularization_budget;
    result.regularization_contribution = round.regularization_contribution;
    result.regularization_anchor_sink_count =
        round.regularization_anchor_sink_count;
    result.regularization_active_sink_count =
        round.regularization_active_sink_count;
    for (const auto &endpoint : package.constraint_endpoints) {
      result.constrained_labels.push_back(
          mcpd3::ConstraintLabel{endpoint.constraint_id,
                                 endpoint.global_node_id,
                                 endpoint.local_index, round.label});
    }
    return result;
  }

  void scaleObjective(long factor,
                      bool saturate_capacity_overflow = false) override {
    saturate_scale_objective_.push_back(saturate_capacity_overflow);
    scale_factors_.push_back(factor);
  }

  const std::vector<mcpd3::PartitionSolveRequest> &requests() const {
    return requests_;
  }

  const std::vector<long> &scaleFactors() const { return scale_factors_; }
  const std::vector<bool> &saturateScaleObjective() const {
    return saturate_scale_objective_;
  }

private:
  const mcpd3::PartitionPackage &packageForRequest(
      const mcpd3::PartitionSolveRequest &request) const {
    if (request.partition_id < 0) {
      require(packages_.size() == 1,
              "scripted worker needs partition id with multiple packages");
      return packages_.begin()->second;
    }
    auto find_iter = packages_.find(request.partition_id);
    require(find_iter != packages_.end(),
            "scripted worker received unknown partition id");
    return find_iter->second;
  }

  std::map<int, mcpd3::PartitionPackage> packages_;
  std::deque<ScriptedRound> script_;
  std::vector<mcpd3::PartitionSolveRequest> requests_;
  std::vector<long> scale_factors_;
  std::vector<bool> saturate_scale_objective_;
};

class ConcurrencyProbeWorker final : public mcpd3::PartitionWorker {
public:
  ConcurrencyProbeWorker(std::atomic<int> *active_solves,
                         std::atomic<int> *max_active_solves)
      : active_solves_(active_solves),
        max_active_solves_(max_active_solves) {}

  void loadPartition(const mcpd3::PartitionPackage &package) override {
    package_ = package;
  }

  mcpd3::PartitionSolveResult solveRound(
      const mcpd3::PartitionSolveRequest &request) override {
    const int active = active_solves_->fetch_add(1) + 1;
    int previous_max = max_active_solves_->load();
    while (active > previous_max &&
           !max_active_solves_->compare_exchange_weak(previous_max, active)) {
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    active_solves_->fetch_sub(1);

    mcpd3::PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = package_.partition_id;
    result.lower_bound = package_.partition_id + 1;
    for (const auto &endpoint : package_.constraint_endpoints) {
      result.constrained_labels.push_back(
          mcpd3::ConstraintLabel{endpoint.constraint_id,
                                 endpoint.global_node_id,
                                 endpoint.local_index, 0});
    }
    return result;
  }

  void scaleObjective(long factor,
                      bool saturate_capacity_overflow = false) override {
    scale_factor_ = factor;
    saturate_scale_objective_ = saturate_capacity_overflow;
  }

private:
  mcpd3::PartitionPackage package_;
  std::atomic<int> *active_solves_;
  std::atomic<int> *max_active_solves_;
  long scale_factor_ = 1;
  bool saturate_scale_objective_ = false;
};

class BatchRecordingWorker final : public mcpd3::PartitionWorker {
public:
  enum class Mode {
    NORMAL,
    DROP_LAST_RESULT,
    RETURN_UNOWNED_PARTITION,
  };

  explicit BatchRecordingWorker(Mode mode = Mode::NORMAL) : mode_(mode) {}

  explicit BatchRecordingWorker(
      mcpd3::PartitionWorkerResourceEstimate resources)
      : resources_(resources) {}

  BatchRecordingWorker(Mode mode,
                       mcpd3::PartitionWorkerResourceEstimate resources)
      : mode_(mode), resources_(resources) {}

  mcpd3::PartitionWorkerResourceEstimate resourceEstimate() const override {
    return resources_;
  }

  void loadPartition(const mcpd3::PartitionPackage &package) override {
    packages_[package.partition_id] = package;
    loaded_partition_ids_.push_back(package.partition_id);
  }

  mcpd3::PartitionSolveResult solveRound(
      const mcpd3::PartitionSolveRequest &request) override {
    ++single_solve_count_;
    return makeResult(request, request.partition_id);
  }

  std::vector<mcpd3::PartitionSolveResult> solveRoundBatch(
      const std::vector<mcpd3::PartitionSolveRequest> &requests) override {
    batch_sizes_.push_back(requests.size());
    std::vector<mcpd3::PartitionSolveResult> results;
    results.reserve(requests.size());
    for (const auto &request : requests) {
      results.push_back(makeResult(request, request.partition_id));
    }
    if (mode_ == Mode::DROP_LAST_RESULT && !results.empty()) {
      results.pop_back();
    }
    if (mode_ == Mode::RETURN_UNOWNED_PARTITION && !results.empty()) {
      results.front().partition_id =
          results.front().partition_id == 0 ? 1 : 0;
    }
    return results;
  }

  void scaleObjective(long, bool = false) override {}

  const std::vector<size_t> &batchSizes() const { return batch_sizes_; }
  int singleSolveCount() const { return single_solve_count_; }
  const std::vector<int> &loadedPartitionIds() const {
    return loaded_partition_ids_;
  }

private:
  mcpd3::PartitionSolveResult makeResult(
      const mcpd3::PartitionSolveRequest &request, int partition_id) const {
    const auto find_iter = packages_.find(request.partition_id);
    require(find_iter != packages_.end(),
            "batch worker received unknown partition id");
    const auto &package = find_iter->second;

    mcpd3::PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = partition_id;
    result.lower_bound = package.partition_id + 1;
    for (const auto &endpoint : package.constraint_endpoints) {
      result.constrained_labels.push_back(
          mcpd3::ConstraintLabel{endpoint.constraint_id,
                                 endpoint.global_node_id,
                                 endpoint.local_index, 0});
    }
    return result;
  }

  Mode mode_ = Mode::NORMAL;
  mcpd3::PartitionWorkerResourceEstimate resources_;
  std::map<int, mcpd3::PartitionPackage> packages_;
  std::vector<int> loaded_partition_ids_;
  std::vector<size_t> batch_sizes_;
  int single_solve_count_ = 0;
};

mcpd3::PartitionPackage makeCoordinatorPackage(int partition_id,
                                               bool is_source) {
  mcpd3::PartitionPackage package;
  package.partition_id = partition_id;
  package.local_node_count = 1;
  package.terminal_capacities = {0};
  package.local_to_global = {11};
  package.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/7,
                                        /*global_node_id=*/11,
                                        /*local_index=*/0,
                                        /*is_source=*/is_source,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});
  return package;
}

std::vector<mcpd3::PartitionPackage> makeCoordinatorPackages() {
  return {makeCoordinatorPackage(/*partition_id=*/0, /*is_source=*/true),
          makeCoordinatorPackage(/*partition_id=*/1, /*is_source=*/false)};
}

void coordinatorCollectsFinalLabelsWhenRequested() {
  const auto scratch =
      std::filesystem::temp_directory_path() /
      ("mcpd3-coordinator-label-test-" +
       std::to_string(std::chrono::steady_clock::now()
                          .time_since_epoch()
                          .count()));
  std::filesystem::create_directories(scratch);
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());

  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 4;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.collect_final_labels = true;
  options.final_label_storage.mode =
      mcpd3::SolverStorageMode::FILE_BACKED_MMAP;
  options.final_label_storage.directory = scratch.string();

  mcpd3::PartitionWorkerCoordinator coordinator(
      makeCoordinatorPackages(), std::move(workers), options);
  const auto result = coordinator.solve();

  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "coordinator should solve the agreeing toy problem");
  require(result.final_labels.size() == 2,
          "coordinator should collect one full label per local node copy");
  require(result.final_labels.isFileBacked(),
          "coordinator final labels should honor mapped result storage");
  require(result.final_labels[0].global_node_id == 11 &&
              result.final_labels[1].global_node_id == 11,
          "coordinator final labels should preserve global node ids");
  require(result.final_labels[0].local_index == 0 &&
              result.final_labels[1].local_index == 0,
          "coordinator final labels should preserve local indices");
  require(result.final_labels[0].label == result.final_labels[1].label,
          "coordinator final labels should agree at convergence");
  std::filesystem::remove_all(scratch);
}

mcpd3::PartitionWorkerResourceEstimate makeWorkerResources(int cpu_count,
                                                           long ram_gb) {
  mcpd3::PartitionWorkerResourceEstimate resources;
  resources.cpu_count = cpu_count;
  resources.ram_gb = ram_gb;
  return resources;
}

mcpd3::PartitionPackage makeSizedPackage(int partition_id,
                                         int local_node_count) {
  mcpd3::PartitionPackage package;
  package.partition_id = partition_id;
  package.local_node_count = local_node_count;
  package.terminal_capacities.assign(static_cast<size_t>(local_node_count), 0);
  package.local_to_global.reserve(static_cast<size_t>(local_node_count));
  for (int i = 0; i < local_node_count; ++i) {
    package.local_to_global.push_back(partition_id * 1000 + i);
  }
  return package;
}

std::vector<mcpd3::PartitionPackage> makeTieBreakRegularizationPackages() {
  mcpd3::PartitionPackage source;
  source.partition_id = 0;
  source.local_node_count = 1;
  source.terminal_capacities = {-10};
  source.local_to_global = {11};
  source.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/7,
                                        /*global_node_id=*/11,
                                        /*local_index=*/0,
                                        /*is_source=*/true,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});

  mcpd3::PartitionPackage target;
  target.partition_id = 1;
  target.local_node_count = 1;
  target.terminal_capacities = {10};
  target.local_to_global = {11};
  target.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/7,
                                        /*global_node_id=*/11,
                                        /*local_index=*/0,
                                        /*is_source=*/false,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});

  return {source, target};
}

std::vector<mcpd3::PartitionPackage> makeOppositeDirectionCyclePackages(
    int source_terminal_capacity = -10, int target_terminal_capacity = 8) {
  mcpd3::PartitionPackage source;
  source.partition_id = 0;
  source.local_node_count = 1;
  source.terminal_capacities = {source_terminal_capacity};
  source.local_to_global = {11};
  source.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/7,
                                        /*global_node_id=*/11,
                                        /*local_index=*/0,
                                        /*is_source=*/true,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});

  mcpd3::PartitionPackage target;
  target.partition_id = 1;
  target.local_node_count = 1;
  target.terminal_capacities = {target_terminal_capacity};
  target.local_to_global = {11};
  target.constraint_endpoints.push_back(
      mcpd3::ConstraintEndpointBinding{/*constraint_id=*/7,
                                        /*global_node_id=*/11,
                                        /*local_index=*/0,
                                        /*is_source=*/false,
                                        /*alpha=*/0,
                                        /*last_alpha=*/0,
                                        /*alpha_momentum=*/0});

  return {source, target};
}

mcpd3::PartitionWorkerCoordinator makeScriptedCoordinator(
    std::deque<ScriptedRound> source_script,
    std::deque<ScriptedRound> target_script,
    const mcpd3::PartitionWorkerCoordinatorOptions &options,
    ScriptedPartitionWorker **source_worker = nullptr,
    ScriptedPartitionWorker **target_worker = nullptr) {
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  auto source = std::make_unique<ScriptedPartitionWorker>(
      std::move(source_script));
  auto target = std::make_unique<ScriptedPartitionWorker>(
      std::move(target_script));
  if (source_worker != nullptr) {
    *source_worker = source.get();
  }
  if (target_worker != nullptr) {
    *target_worker = target.get();
  }
  workers.push_back(std::move(source));
  workers.push_back(std::move(target));
  return mcpd3::PartitionWorkerCoordinator(makeCoordinatorPackages(),
                                           std::move(workers), options);
}

void coordinatorScaledEpsilonControlsMatchNativeSchedule() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 12;
  options.max_iteration_count = 2;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 100;
  options.scaled_epsilon_max_step_size = 12;
  options.scaled_epsilon_strength_cap = 3;

  ScriptedPartitionWorker *source_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 0}}, options,
      &source_worker);
  const auto result = coordinator.solve();

  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "custom scaled-epsilon coordinator schedule should agree");
  require(source_worker->requests().size() == 2,
          "custom scaled-epsilon schedule should run two rounds");
  for (const auto &request : source_worker->requests()) {
    require(request.regularization_strength == 3,
            "coordinator must apply the configured scaled-epsilon cap");
  }

  options.scaled_epsilon_max_step_size = 11;
  ScriptedPartitionWorker *disabled_source_worker = nullptr;
  auto disabled = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 0}}, options,
      &disabled_source_worker);
  (void)disabled.solve();
  for (const auto &request : disabled_source_worker->requests()) {
    require(request.regularization_strength == 0,
            "coordinator must honor the configured scaled-epsilon cutoff");
  }
}

void coordinatorDisagreementPlateauMatchesNativePulses() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 6;
  options.max_total_iteration_count = 6;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.disagreement_patience = 2;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 1000;
  options.regularization_scheme =
      mcpd3::PartitionWorkerRegularizationScheme::
          DISAGREEMENT_PLATEAU_EPSILON;

  ScriptedPartitionWorker *source_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}, {10, 0}, {10, 0},
                                {10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}, {20, 1}, {20, 1},
                                {20, 1}, {20, 0}},
      options, &source_worker);
  const auto result = coordinator.solve();

  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "plateau regularization fixture should eventually agree");
  require(result.disagreement_plateau_activation_count == 1,
          "persistent disagreement should activate plateau mode once");
  const std::vector<mcpd3::Capacity> expected_strengths{0, 0, 0, 1, 0, 1};
  require(source_worker->requests().size() == expected_strengths.size(),
          "plateau fixture request count mismatch");
  for (size_t i = 0; i < expected_strengths.size(); ++i) {
    require(source_worker->requests()[i].regularization_strength ==
                expected_strengths[i],
            "plateau regularization must pulse only after each complete flat "
            "window");
  }
}

void coordinatorHonorsTotalIterationBudget() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 10;
  options.max_total_iteration_count = 2;
  options.num_optimization_scales = 3;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.promote_objective_scale_on_overbudget = false;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}}, options);
  const auto result = coordinator.solve();
  require(result.total_iterations == 2,
          "coordinator must stop at the global iteration budget");
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  ITERATION_COUNT_EXCEEDED,
          "global iteration exhaustion must not report optimality");
}

void coordinatorCanRetryUnitStepWithoutMomentum() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 1;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = true;
  options.retry_unit_step_without_momentum = true;
  options.promote_objective_scale_on_overbudget = false;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 0}}, options);
  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "unit-step no-momentum retry should reach agreement");
  require(result.unit_step_no_momentum_retry_count == 1,
          "coordinator should perform exactly one unit-step cleanup retry");
  require(result.total_iterations == 2,
          "unit-step cleanup retry should add one solve round");
}

void coordinatorRejectsInvalidNativeScheduleControls() {
  auto require_invalid = [](const mcpd3::PartitionWorkerCoordinatorOptions &opts,
                            const std::string &message) {
    requireThrows(
        [&] {
          auto coordinator = makeScriptedCoordinator(
              std::deque<ScriptedRound>{{10, 0}},
              std::deque<ScriptedRound>{{20, 0}}, opts);
          (void)coordinator;
        },
        message);
  };

  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.max_total_iteration_count = -1;
  require_invalid(options, "negative global iteration budget must be rejected");
  options.max_total_iteration_count = 0;
  options.scaled_epsilon_max_step_size = 0;
  require_invalid(options, "zero scaled-epsilon cutoff must be rejected");
  options.scaled_epsilon_max_step_size = 10;
  options.scaled_epsilon_strength_cap = -1;
  require_invalid(options, "negative scaled-epsilon cap must be rejected");
  options.scaled_epsilon_strength_cap = 0;
  options.regularization_scheme =
      mcpd3::PartitionWorkerRegularizationScheme::
          DISAGREEMENT_PLATEAU_EPSILON;
  options.disagreement_patience = 0;
  require_invalid(options, "zero plateau patience must be rejected");
}

void coordinatorRunRoundReportsCertifiedRegularizedLowerBound() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.objective_scale = 100;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{100, 0, 10, 7, 1, 1}},
      std::deque<ScriptedRound>{{50, 1, 5, 2, 1, 1}}, options);

  const auto stats = coordinator.runRound(
      /*round_id=*/1, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/10);

  require(stats.original_objective == 150,
          "original objective should sum solver-reported original values");
  require(stats.regularized_objective == 159,
          "regularized objective should be original value plus contribution");
  require(stats.certified_lower_bound == 144,
          "certified field should subtract full budget from regularized objective");
  require(stats.lower_bound == 144,
          "certified lower bound should subtract full budget from "
          "regularized objective");
  require(stats.regularization_budget == 15,
          "round should sum regularization budgets");
  require(stats.regularization_contribution == 9,
          "round should sum actual regularization contributions");
  require(stats.disagreement_count == 1,
          "generic regularized certificate test should remain infeasible");
}

void coordinatorRunRoundStrengthensCertificateAtAgreement() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.objective_scale = 100;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{100, 0, 10, 7, 1, 1}},
      std::deque<ScriptedRound>{{50, 0, 5, 2, 1, 1}}, options);

  const auto stats = coordinator.runRound(
      /*round_id=*/1, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/10);

  require(stats.disagreement_count == 0,
          "agreement certificate test should be globally feasible");
  require(stats.regularization_budget == 15 &&
              stats.regularization_budget < options.objective_scale,
          "agreement certificate test requires strict secondary budget");
  require(stats.certified_lower_bound == stats.original_objective,
          "agreement under a strict lexicographic budget should certify the "
          "primary objective");
  require(stats.lower_bound == stats.original_objective,
          "selected lower bound should use the strengthened agreement "
          "certificate");
}

void coordinatorRunRoundLeavesUnregularizedLowerBoundUnchanged() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{12, 0}},
      std::deque<ScriptedRound>{{30, 0}}, options);

  const auto stats = coordinator.runRound(
      /*round_id=*/1, /*scale=*/100, /*step_size=*/100,
      /*regularization_strength=*/0);

  require(stats.original_objective == 42,
          "unregularized original objective should equal solver value");
  require(stats.certified_lower_bound == 42,
          "unregularized certified field should equal selected value");
  require(stats.lower_bound == 42,
          "unregularized certified lower bound should equal selected value");
  require(stats.regularized_objective == 42,
          "unregularized regularized objective diagnostic should match");
  require(stats.regularization_budget == 0,
          "unregularized round should have zero budget");
  require(stats.regularization_contribution == 0,
          "unregularized round should have zero contribution");
}

std::vector<std::unique_ptr<mcpd3::PartitionWorker>> makeInProcessWorkers(
    size_t count) {
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  for (size_t i = 0; i < count; ++i) {
    workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  }
  return workers;
}

struct OneNodeSourceSolveResult {
  mcpd3::Objective lower_bound = 0;
  int label = 0;
  mcpd3::Objective regularization_budget = 0;
  mcpd3::Objective regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
};

OneNodeSourceSolveResult solveOneNodeSourceProblem(
    long alpha, long last_alpha, int terminal_capacity,
    int regularization_strength, bool seed_sink_anchor) {
  std::list<mcpd3::DualDecompositionConstraintArc> constraints;
  constraints.emplace_back(alpha, last_alpha, /*alpha_momentum=*/0,
                           /*partition_index_source=*/0,
                           /*partition_index_target=*/1,
                           /*local_index_source=*/0,
                           /*local_index_target=*/-1);
  auto ref = --constraints.end();
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{terminal_capacity});
  solver.addSourceDualDecompositionConstraint(ref);
  if (seed_sink_anchor) {
    solver.setMinCutSolution(std::vector<int>{1});
  }
  solver.setRegularizationStrength(regularization_strength);
  solver.solve();
  return OneNodeSourceSolveResult{
      solver.getMinCutValue(),
      solver.getMinCutSolution(/*index=*/0),
      solver.getLastRegularizationBudget(),
      solver.getLastRegularizationContribution(),
      solver.getLastRegularizationAnchorSinkCount(),
      solver.getLastRegularizationActiveSinkCount()};
}

void scaledEpsilonRegularizationUsesPreviousSinkAnchors() {
  const auto no_anchor =
      solveOneNodeSourceProblem(/*alpha=*/-100, /*last_alpha=*/-90,
                                /*terminal_capacity=*/-100,
                                /*regularization_strength=*/10,
                                /*seed_sink_anchor=*/false);
  require(no_anchor.regularization_budget == 0,
          "regularization should not apply without an existing anchor");
  require(no_anchor.lower_bound == 100,
          "unanchored tie lower bound should be unregularized");

  const auto unchanged_alpha =
      solveOneNodeSourceProblem(/*alpha=*/-100, /*last_alpha=*/-100,
                                /*terminal_capacity=*/-100,
                                /*regularization_strength=*/10,
                                /*seed_sink_anchor=*/true);
  require(unchanged_alpha.regularization_budget == 0,
          "unchanged alpha terms should not activate sink anchors");
  require(unchanged_alpha.lower_bound == 100,
          "unchanged-alpha lower bound should remain unregularized");

  const auto strict =
      solveOneNodeSourceProblem(/*alpha=*/-89, /*last_alpha=*/-100,
                                /*terminal_capacity=*/-100,
                                /*regularization_strength=*/10,
                                /*seed_sink_anchor=*/true);
  require(strict.label == 1,
          "scaled epsilon regularization must preserve larger strict gaps");
  require(strict.lower_bound == 89,
          "strict optimum lower bound should exclude regularization");
  require(strict.regularization_budget == 10,
          "strict anchored solve should report the scaled epsilon budget");
  require(strict.regularization_contribution == 10,
          "strict anchored sink solution should pay the regularizer");
  require(strict.regularization_anchor_sink_count == 1,
          "strict anchored solve should count the sink anchor");
  require(strict.regularization_active_sink_count == 1,
          "strict anchored sink solution should count active regularization");

  const auto tied =
      solveOneNodeSourceProblem(/*alpha=*/-100, /*last_alpha=*/-90,
                                /*terminal_capacity=*/-100,
                                /*regularization_strength=*/10,
                                /*seed_sink_anchor=*/true);
  require(tied.label == 0,
          "scaled epsilon regularization should select source on ties");
  require(tied.lower_bound == 100,
          "tie-broken lower bound should remain the unregularized optimum");
  require(tied.regularization_budget == 10,
          "tie-broken solve should report the scaled epsilon budget");
  require(tied.regularization_contribution == 0,
          "tie-broken source solution should not pay regularization");
  require(tied.regularization_anchor_sink_count == 1,
          "tie-broken solve should count the sink anchor");
  require(tied.regularization_active_sink_count == 0,
          "tie-broken source solution should not count active regularization");
}

void scaledEpsilonRegularizationPersistsAfterSourceTieBreak() {
  std::list<mcpd3::DualDecompositionConstraintArc> constraints;
  constraints.emplace_back(/*alpha=*/-100, /*last_alpha=*/-90,
                           /*alpha_momentum=*/0,
                           /*partition_index_source=*/0,
                           /*partition_index_target=*/1,
                           /*local_index_source=*/0,
                           /*local_index_target=*/-1);
  auto ref = --constraints.end();
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{-100});
  solver.addSourceDualDecompositionConstraint(ref);
  solver.setMinCutSolution(std::vector<int>{1});
  solver.setRegularizationStrength(10);

  solver.solve();
  require(solver.getMinCutSolution(/*index=*/0) == 0,
          "changed-alpha epsilon should break the initial tie to source");
  require(solver.getLastRegularizationBudget() == 10,
          "changed-alpha solve should activate the epsilon budget");
  require(solver.getLastRegularizationContribution() == 0,
          "source label should not pay the active epsilon term");

  ref->last_alpha = ref->alpha;
  solver.solve();
  require(solver.getLastRegularizationBudget() == 10,
          "unchanged-alpha solve should keep the active epsilon term");
  require(solver.getMinCutSolution(/*index=*/0) == 0,
          "persistent epsilon should keep the tied source label");

  ref->last_alpha = ref->alpha;
  ref->alpha = -90;
  solver.solve();
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
  require(solver.getLastRegularizationBudget() == 0,
          "legacy replay should clear its binary anchor after a source update");
  require(solver.getMinCutSolution(/*index=*/0) == 1,
          "legacy replay should restore its historical unregularized label");
#else
  require(solver.getLastRegularizationBudget() == 10,
          "cumulative epsilon must persist after a source tie-break");
  require(solver.getMinCutSolution(/*index=*/0) == 0,
          "persistent cumulative epsilon should retain the source tie-break");
#endif
}

void scaledEpsilonRegularizationAccumulatesAcrossPersistentSinkUpdates() {
  std::list<mcpd3::DualDecompositionConstraintArc> constraints;
  constraints.emplace_back(/*alpha=*/-989, /*last_alpha=*/-1000,
                           /*alpha_momentum=*/0,
                           /*partition_index_source=*/0,
                           /*partition_index_target=*/1,
                           /*local_index_source=*/0,
                           /*local_index_target=*/-1);
  auto ref = --constraints.end();
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{-1000});
  solver.addSourceDualDecompositionConstraint(ref);
  solver.setMinCutSolution(std::vector<int>{1});
  solver.setRegularizationStrength(10);

  solver.solve();
  require(solver.getMinCutSolution(/*index=*/0) == 1,
          "first cumulative epsilon update should preserve a strict sink");
  require(solver.getLastRegularizationBudget() == 10,
          "first sink update should consume one epsilon unit");

  ref->last_alpha = ref->alpha;
  ref->alpha = -978;
  solver.solve();
  require(solver.getMinCutSolution(/*index=*/0) == 1,
          "second cumulative epsilon update should preserve a strict sink");
  require(solver.getLastRegularizationBudget() == 20,
          "persistent sink disagreement should accumulate epsilon budget");

  ref->last_alpha = ref->alpha;
  solver.solve();
  require(solver.getLastRegularizationBudget() == 20,
          "unchanged alpha must not consume additional epsilon budget");

  solver.setRegularizationStrength(1);
  ref->last_alpha = ref->alpha;
  ref->alpha = -967;
  solver.solve();
  require(solver.getMinCutSolution(/*index=*/0) == 1,
          "unit-scale cumulative update should preserve a strict sink");
  require(solver.getLastRegularizationBudget() == 21,
          "unit-scale epsilon must add to the prior high-scale budget");
  require(solver.getLastRegularizationContribution() == 21,
          "active sink should pay the complete cumulative regularizer");
}

void lowScaleScaledEpsilonRegularizationHandlesBoundaryTie() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.min_step_size = 1;
  options.objective_scale = 10000;

  const auto packages = makeTieBreakRegularizationPackages();
  mcpd3::PartitionWorkerCoordinator unregularized_coordinator(
      packages, makeInProcessWorkers(packages.size()), options);
  const auto initial_unregularized = unregularized_coordinator.runRound(
      /*round_id=*/1, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/0);
  require(initial_unregularized.disagreement_count == 1,
          "initial unregularized round should disagree");
  require(initial_unregularized.regularization_budget == 0,
          "unregularized round should not report regularization budget");

  const auto tied_unregularized = unregularized_coordinator.runRound(
      /*round_id=*/2, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/0);
  require(tied_unregularized.disagreement_count == 0,
          "current unregularized maxflow tie-break should choose source");
  require(tied_unregularized.regularization_budget == 0,
          "unregularized tie should not report regularization budget");

  mcpd3::PartitionWorkerCoordinator regularized_coordinator(
      packages, makeInProcessWorkers(packages.size()), options);
  const auto regularized_initial = regularized_coordinator.runRound(
      /*round_id=*/1, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/10);
  require(regularized_initial.disagreement_count == 1,
          "first low-scale round should still disagree before anchors exist");
  require(regularized_initial.regularization_budget == 0,
          "first low-scale round should not regularize without anchors");

  const auto regularized = regularized_coordinator.runRound(
      /*round_id=*/2, /*scale=*/10, /*step_size=*/10,
      /*regularization_strength=*/10);
  require(regularized.disagreement_count == 0,
          "scaled epsilon regularization should select agreeing tied labels");
  require(regularized.regularization_budget == 10,
          "regularized round should anchor the source-side sink label");
  require(regularized.regularization_contribution == 0,
          "agreeing solution should not pay regularization");
}

void fullSolveStopsOptimalOnUnregularizedAgreement() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}},
      std::deque<ScriptedRound>{{15, 0}}, options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "unregularized agreement should be optimal");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::NO_DISAGREEMENT,
          "optimal agreement should report no-disagreement stop");
  require(result.best_lower_bound_raw == 25,
          "best lower bound should include all partitions");
  require(result.final_disagreement_count == 0,
          "final disagreement count should be zero");
  require(result.total_iterations == 1, "optimal solve should stop in one iter");
  require(result.progress_records.size() == 1,
          "progress should record every iteration");
  const auto &record = result.progress_records.front();
  require(record.scale == 100, "progress record scale mismatch");
  require(record.iteration == 0, "progress record iteration mismatch");
  require(record.total_iteration == 1,
          "progress record total iteration mismatch");
  require(record.max_iteration == 5, "progress record max iteration mismatch");
  require(record.lower_bound == 25, "progress record lower bound mismatch");
  require(record.best_lower_bound == 25,
          "progress record best lower bound mismatch");
  require(record.certified_lower_bound == 25,
          "progress record certified lower bound mismatch");
  require(record.best_certified_lower_bound == 25,
          "progress record best certified lower bound mismatch");
  require(record.regularized_objective == 25,
          "unregularized progress diagnostic should match lower bound");
  require(record.best_regularized_objective == 25,
          "unregularized best regularized diagnostic should match lower bound");
  require(record.regularization_strength == 0,
          "progress record regularization strength mismatch");
}

void fullSolveReportsProgressThroughConfiguredCallback() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.progress_report_interval = 2;

  std::vector<mcpd3::PartitionWorkerProgressRecord> callbacks;
  options.progress_callback =
      [&](const mcpd3::PartitionWorkerProgressRecord &record) {
        callbacks.push_back(record);
      };

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}, {12, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 1}, {22, 0}}, options);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "callback progress solve should finish on agreement");
  require(result.total_iterations == 3,
          "callback progress solve should run three iterations");
  require(result.progress_records.size() == 3,
          "callback progress should not replace stored progress records");
  require(callbacks.size() == 1,
          "callback should fire only on configured interval");
  require(callbacks[0].total_iteration == 2,
          "callback should report the second total iteration");
  require(callbacks[0].lower_bound == 32,
          "callback should include the round lower bound");
  require(callbacks[0].best_lower_bound == 32,
          "callback should include the best lower bound");
  require(callbacks[0].certified_lower_bound == 32,
          "callback should include the certified lower bound");
  require(callbacks[0].best_certified_lower_bound == 32,
          "callback should include the best certified lower bound");
  require(callbacks[0].regularized_objective == 32,
          "callback should include the regularized objective diagnostic");
  require(callbacks[0].best_regularized_objective == 32,
          "callback should include the best regularized objective diagnostic");
  require(callbacks[0].disagreement_count == 1,
          "callback should include disagreement count");
}

void disabledProgressCallbackDoesNotFire() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.progress_report_interval = 0;

  bool called = false;
  options.progress_callback =
      [&](const mcpd3::PartitionWorkerProgressRecord &) { called = true; };

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}},
      std::deque<ScriptedRound>{{20, 0}}, options);
  (void)coordinator.solve();
  require(!called, "disabled progress interval should not call callback");
}

void invalidProgressIntervalIsRejected() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.progress_report_interval = -1;
  bool threw = false;
  try {
    auto coordinator = makeScriptedCoordinator(
        std::deque<ScriptedRound>{{10, 0}},
        std::deque<ScriptedRound>{{20, 0}}, options);
  } catch (const std::runtime_error &e) {
    threw = std::string(e.what()).find("progress report interval") !=
            std::string::npos;
  }
  require(threw, "negative progress interval should be rejected");
}

void fullSolveReportsBestBoundUsingObjectiveScale() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.objective_scale = 1;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}},
      std::deque<ScriptedRound>{{15, 0}}, options);

  const auto result = coordinator.solve();
  require(result.best_lower_bound_raw == 25,
          "best lower bound raw value should include all partitions");
  require(result.best_lower_bound == 25,
          "reported best lower bound should not be divided by step size");

  options.objective_scale = 5;
  auto scaled_coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{50, 0}},
      std::deque<ScriptedRound>{{75, 0}}, options);
  const auto scaled_result = scaled_coordinator.solve();
  require(scaled_result.best_lower_bound_raw == 125,
          "scaled raw bound should include all partitions");
  require(scaled_result.best_lower_bound == 25,
          "explicit objective scale should control reported lower bound");

  options.objective_scale = 0;
  bool threw = false;
  try {
    auto invalid_coordinator = makeScriptedCoordinator(
        std::deque<ScriptedRound>{{10, 0}},
        std::deque<ScriptedRound>{{15, 0}}, options);
  } catch (const std::runtime_error &e) {
    threw = std::string(e.what()).find("objective scale") !=
            std::string::npos;
  }
  require(threw, "nonpositive coordinator objective scale should fail");
}

void fullSolveReportsLowScaleRegularizedAgreementAsOptimal() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10000;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0, 2, 1, 3, 1}},
      std::deque<ScriptedRound>{{15, 0, 4, 2, 5, 2}}, options,
      &source_worker, &target_worker);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "scaled epsilon regularized agreement should certify optimality");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "regularized agreement should report its stop reason");
  require(result.total_iterations == 1,
          "solve should not run an unregularized confirmation round");
  require(result.final_objective_raw == 25,
          "regularized agreement should preserve the final objective");
  require(result.final_certified_lower_bound_raw == 25,
          "regularized agreement should certify the feasible primary objective");
  require(result.final_regularized_objective_raw == 28,
          "regularized agreement should preserve the final regularized objective");
  require(result.best_lower_bound_raw == 25,
          "regularized agreement should store the exact certified lower bound");
  require(result.best_certified_lower_bound_raw == 25,
          "regularized agreement should store the explicit certified lower bound");
  require(result.best_regularized_objective_raw == 28,
          "regularized agreement should preserve the regularized objective");
  require(result.final_objective > 0.002499 &&
              result.final_objective < 0.002501,
          "regularized final objective should use objective scale");
  require(result.best_lower_bound > 0.002499 &&
              result.best_lower_bound < 0.002501,
          "regularized certified lower bound should use objective scale");
  require(result.best_certified_lower_bound > 0.002499 &&
              result.best_certified_lower_bound < 0.002501,
          "explicit certified lower bound should use objective scale");
  require(result.best_regularized_objective > 0.002799 &&
              result.best_regularized_objective < 0.002801,
          "regularized objective diagnostic should use objective scale");
  require(result.progress_records.size() == 1,
          "regularized agreement should record one progress row");
  require(result.progress_records[0].lower_bound == 25,
          "regularized progress should report certified lower bound");
  require(result.progress_records[0].best_lower_bound == 25,
          "regularized progress should report best certified lower bound");
  require(result.progress_records[0].certified_lower_bound == 25,
          "regularized progress should report certified lower bound explicitly");
  require(result.progress_records[0].best_certified_lower_bound == 25,
          "regularized progress should report best certified lower bound explicitly");
  require(result.progress_records[0].regularized_objective == 28,
          "regularized progress should report original plus contribution");
  require(result.progress_records[0].best_regularized_objective == 28,
          "regularized progress should report best regularized objective");
  require(result.progress_records[0].regularization_strength == 10,
          "regularized round should use regularization");
  require(result.progress_records[0].regularization_budget == 6,
          "regularized round should report regularization budget");
  require(result.progress_records[0].disagreement_count == 0,
          "regularized round should reach agreement");
  require(result.final_regularization_budget == 6,
          "final diagnostics should come from the scaled epsilon round");
  require(source_worker->requests()[0].regularization_strength == 10,
          "first source request should be regularized");
  require(source_worker->requests().size() == 1,
          "source worker should not receive a confirmation request");
  require(target_worker->requests().size() == 1,
          "target worker should not receive a confirmation request");
}

void fullSolveUsesScaledEpsilonRegularizationToReachAgreement() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10000;

  const auto packages = makeTieBreakRegularizationPackages();
  mcpd3::PartitionWorkerCoordinator coordinator(
      packages, makeInProcessWorkers(packages.size()), options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "scaled epsilon tie-break agreement should be optimal");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "scaled epsilon tie-break stop reason mismatch");
  require(result.total_iterations == 2,
          "tie-break solve should stop after the anchored low-scale round");
  require(result.final_regularization_budget == 10,
          "final diagnostics should include the tie-break budget");
  require(result.final_regularization_contribution == 0,
          "agreeing tie-break should not pay the regularizer");
  require(result.progress_records.size() == 2,
          "tie-break solve should record initial and regularized rounds");
  require(result.progress_records[0].regularization_strength == 10,
          "low-scale schedule should request regularization immediately");
  require(result.progress_records[0].regularization_budget == 0,
          "first low-scale round should have no anchors yet");
  require(result.progress_records[0].disagreement_count == 1,
          "first round should expose the initial disagreement");
  require(result.progress_records[1].regularization_strength == 10,
          "second low-scale round should remain regularized");
  require(result.progress_records[1].regularization_budget == 10,
          "second low-scale round should report its regularization budget");
  require(result.progress_records[1].disagreement_count == 0,
          "second low-scale round should reach agreement");
}

void coordinatorRoutesMultiplePackagesToOneWorker() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 2;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;

  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  auto scripted_worker = std::make_unique<ScriptedPartitionWorker>(
      std::deque<ScriptedRound>{{10, 0}, {20, 1}, {30, 1}, {40, 1}});
  auto *worker_ptr = scripted_worker.get();
  workers.push_back(std::move(scripted_worker));

  mcpd3::PartitionWorkerCoordinator coordinator(
      makeCoordinatorPackages(), std::move(workers), options);

  const auto first =
      coordinator.runRound(/*round_id=*/1, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);
  const auto second =
      coordinator.runRound(/*round_id=*/2, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);

  require(first.lower_bound == 30,
          "single worker should solve all first-round packages");
  require(first.disagreement_count == 1,
          "first single-worker round should see disagreement");
  require(second.lower_bound == 70,
          "single worker should solve all second-round packages");
  require(second.disagreement_count == 0,
          "second single-worker round should reach agreement");
  require(worker_ptr->requests().size() == 4,
          "single worker should receive one request per package per round");
  require(worker_ptr->requests()[0].partition_id == 0,
          "first request should target partition 0");
  require(worker_ptr->requests()[1].partition_id == 1,
          "second request should target partition 1");
  require(worker_ptr->requests()[2].partition_id == 0,
          "third request should return to partition 0");
  require(worker_ptr->requests()[3].partition_id == 1,
          "fourth request should return to partition 1");
  require(worker_ptr->requests()[2].alpha_updates.size() == 1,
          "partition 0 should receive its alpha update");
  require(worker_ptr->requests()[3].alpha_updates.size() == 1,
          "partition 1 should receive its alpha update");
  require(worker_ptr->requests()[2].alpha_updates[0].alpha == 100,
          "partition 0 alpha should update after disagreement");
  require(worker_ptr->requests()[3].alpha_updates[0].alpha == 100,
          "partition 1 alpha should update after disagreement");
}

void coordinatorBatchesMultiplePackagesPerWorker() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;

  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  auto worker = std::make_unique<BatchRecordingWorker>();
  auto *worker_ptr = worker.get();
  workers.push_back(std::move(worker));

  mcpd3::PartitionWorkerCoordinator coordinator(
      makeCoordinatorPackages(), std::move(workers), options);
  const auto stats =
      coordinator.runRound(/*round_id=*/7, /*scale=*/100,
                           /*step_size=*/100, /*regularization_strength=*/0);

  require(stats.lower_bound == 3,
          "batched worker lower bound should combine both partitions");
  require(stats.disagreement_count == 0,
          "batched worker should return agreeing labels");
  require(worker_ptr->batchSizes() == std::vector<size_t>{2},
          "coordinator should send one batch containing both partitions");
  require(worker_ptr->singleSolveCount() == 0,
          "coordinator should use solveRoundBatch instead of per-partition RPCs");
}

void coordinatorBalancesInitialPackagesByCpuCapacity() {
  std::vector<mcpd3::PartitionPackage> packages{
      makeSizedPackage(/*partition_id=*/0, /*local_node_count=*/100),
      makeSizedPackage(/*partition_id=*/1, /*local_node_count=*/10),
      makeSizedPackage(/*partition_id=*/2, /*local_node_count=*/10),
      makeSizedPackage(/*partition_id=*/3, /*local_node_count=*/10),
      makeSizedPackage(/*partition_id=*/4, /*local_node_count=*/10),
      makeSizedPackage(/*partition_id=*/5, /*local_node_count=*/10)};

  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  auto small_worker = std::make_unique<BatchRecordingWorker>(
      makeWorkerResources(/*cpu_count=*/1, /*ram_gb=*/4));
  auto large_worker = std::make_unique<BatchRecordingWorker>(
      makeWorkerResources(/*cpu_count=*/4, /*ram_gb=*/16));
  auto *small_worker_ptr = small_worker.get();
  auto *large_worker_ptr = large_worker.get();
  workers.push_back(std::move(small_worker));
  workers.push_back(std::move(large_worker));

  mcpd3::PartitionWorkerCoordinator coordinator(
      std::move(packages), std::move(workers));

  require(small_worker_ptr->loadedPartitionIds() ==
              std::vector<int>({1, 2}),
          "low-CPU worker should receive only the first small packages");
  require(large_worker_ptr->loadedPartitionIds() ==
              std::vector<int>({0, 3, 4, 5}),
          "high-CPU worker should receive the largest package and more work");
}

void coordinatorUsesRamAsInitialAssignmentTieBreaker() {
  std::vector<mcpd3::PartitionPackage> packages{
      makeSizedPackage(/*partition_id=*/0, /*local_node_count=*/100),
      makeSizedPackage(/*partition_id=*/1, /*local_node_count=*/10),
      makeSizedPackage(/*partition_id=*/2, /*local_node_count=*/10)};

  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  auto low_ram_worker = std::make_unique<BatchRecordingWorker>(
      makeWorkerResources(/*cpu_count=*/2, /*ram_gb=*/4));
  auto high_ram_worker = std::make_unique<BatchRecordingWorker>(
      makeWorkerResources(/*cpu_count=*/2, /*ram_gb=*/64));
  auto *low_ram_worker_ptr = low_ram_worker.get();
  auto *high_ram_worker_ptr = high_ram_worker.get();
  workers.push_back(std::move(low_ram_worker));
  workers.push_back(std::move(high_ram_worker));

  mcpd3::PartitionWorkerCoordinator coordinator(
      std::move(packages), std::move(workers));

  require(low_ram_worker_ptr->loadedPartitionIds() ==
              std::vector<int>({1, 2}),
          "lower-RAM equal-CPU worker should receive the small packages");
  require(high_ram_worker_ptr->loadedPartitionIds() ==
              std::vector<int>({0}),
          "higher-RAM equal-CPU worker should receive the largest package");
}

void coordinatorRejectsMalformedBatchResponses() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;

  {
    std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
    workers.push_back(std::make_unique<BatchRecordingWorker>(
        BatchRecordingWorker::Mode::DROP_LAST_RESULT));
    mcpd3::PartitionWorkerCoordinator coordinator(
        makeCoordinatorPackages(), std::move(workers), options);
    bool threw = false;
    try {
      (void)coordinator.runRound(/*round_id=*/1, /*scale=*/100,
                                 /*step_size=*/100,
                                 /*regularization_strength=*/0);
    } catch (const std::runtime_error &e) {
      threw = std::string(e.what()).find("result count") !=
              std::string::npos;
    }
    require(threw, "coordinator should reject short batch responses");
  }

  {
    std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
    workers.push_back(std::make_unique<BatchRecordingWorker>(
        BatchRecordingWorker::Mode::RETURN_UNOWNED_PARTITION));
    workers.push_back(std::make_unique<BatchRecordingWorker>());
    mcpd3::PartitionWorkerCoordinator coordinator(
        makeCoordinatorPackages(), std::move(workers), options);
    bool threw = false;
    try {
      (void)coordinator.runRound(/*round_id=*/1, /*scale=*/100,
                                 /*step_size=*/100,
                                 /*regularization_strength=*/0);
    } catch (const std::runtime_error &e) {
      threw = std::string(e.what()).find("unowned partition") !=
              std::string::npos;
    }
    require(threw, "coordinator should reject unowned batch results");
  }
}

void inProcessWorkerBatchSolvesDistinctLoadedPartitions() {
  mcpd3::InProcessPartitionWorker worker;
  for (const auto &package : makeCoordinatorPackages()) {
    worker.loadPartition(package);
  }

  mcpd3::PartitionSolveRequest first_request;
  first_request.round_id = 9;
  first_request.partition_id = 0;
  first_request.scale = 100;
  first_request.regularization_strength = 0;

  mcpd3::PartitionSolveRequest second_request = first_request;
  second_request.partition_id = 1;

  const auto results =
      worker.solveRoundBatch({first_request, second_request});
  require(results.size() == 2,
          "in-process batch should return one result per request");
  require(results[0].round_id == 9 && results[1].round_id == 9,
          "in-process batch should preserve round ids");
  require(results[0].partition_id == 0 && results[1].partition_id == 1,
          "in-process batch should preserve request order and partition ids");
  require(results[0].constrained_labels.size() == 1 &&
              results[1].constrained_labels.size() == 1,
          "in-process batch should return constrained labels");

  bool threw = false;
  try {
    (void)worker.solveRoundBatch({first_request, first_request});
  } catch (const std::runtime_error &e) {
    threw = std::string(e.what()).find("distinct partitions") !=
            std::string::npos;
  }
  require(threw, "in-process batch should reject duplicate partitions");
}

void coordinatorSendsOnlyDirtyAlphaUpdates() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 4;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}, {12, 0}, {13, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 0}, {22, 0}, {23, 0}},
      options, &source_worker, &target_worker);

  const auto first =
      coordinator.runRound(/*round_id=*/1, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);
  const auto second =
      coordinator.runRound(/*round_id=*/2, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);
  const auto third =
      coordinator.runRound(/*round_id=*/3, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);
  const auto fourth =
      coordinator.runRound(/*round_id=*/4, /*scale=*/100, /*step_size=*/100,
                           /*regularization_strength=*/0);

  require(first.disagreement_count == 1,
          "first dirty-alpha round should create disagreement");
  require(second.disagreement_count == 0,
          "second dirty-alpha round should agree");
  require(third.disagreement_count == 0,
          "third dirty-alpha round should agree");
  require(fourth.disagreement_count == 0,
          "fourth dirty-alpha round should agree");
  require(source_worker->requests().size() == 4,
          "source worker should receive four requests");
  require(target_worker->requests().size() == 4,
          "target worker should receive four requests");
  require(source_worker->requests()[0].alpha_updates.empty(),
          "initial zero alpha state should not be resent");
  require(target_worker->requests()[0].alpha_updates.empty(),
          "target initial zero alpha state should not be resent");
  require(source_worker->requests()[1].alpha_updates.size() == 1,
          "changed alpha should be sent on the next source request");
  require(target_worker->requests()[1].alpha_updates.size() == 1,
          "changed alpha should be sent on the next target request");
  require(source_worker->requests()[1].alpha_updates[0].alpha == 100,
          "dirty alpha update should carry the changed alpha");
  require(source_worker->requests()[1].alpha_updates[0].last_alpha == 0,
          "dirty alpha update should carry the previous alpha");
  require(source_worker->requests()[2].alpha_updates.size() == 1,
          "last-alpha catch-up should be sent once after agreement");
  require(source_worker->requests()[2].alpha_updates[0].alpha == 100,
          "catch-up alpha should preserve the alpha value");
  require(source_worker->requests()[2].alpha_updates[0].last_alpha == 100,
          "catch-up update should synchronize last_alpha");
  require(target_worker->requests()[2].alpha_updates.size() == 1,
          "target catch-up update should be sent once");
  require(source_worker->requests()[3].alpha_updates.empty(),
          "unchanged synchronized alpha should not be resent");
  require(target_worker->requests()[3].alpha_updates.empty(),
          "unchanged synchronized target alpha should not be resent");
}

void inProcessCoordinatorSolvesWithOneWorkerOwningAllPackages() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10000;

  const auto packages = makeTieBreakRegularizationPackages();
  mcpd3::PartitionWorkerCoordinator coordinator(
      packages, makeInProcessWorkers(/*count=*/1), options);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "single in-process worker should solve all owned packages");
  require(result.final_disagreement_count == 0,
          "single in-process worker solve should finish with agreement");
  require(result.final_regularization_budget == 10,
          "single in-process worker should preserve regularization diagnostics");
}

void unitScaleResolvesOppositeDirectionCycle() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.min_step_size = 1;
  options.objective_scale = 100;

  auto make_coordinator = [&]() {
    auto packages = makeOppositeDirectionCyclePackages(
        /*source_terminal_capacity=*/-100, /*target_terminal_capacity=*/92);
    return mcpd3::PartitionWorkerCoordinator(
        packages, makeInProcessWorkers(packages.size()), options);
  };

  auto scale_ten = make_coordinator();
  mcpd3::Objective previous_regularization_budget = 0;
  for (int round = 1; round <= 4; ++round) {
    const auto stats = scale_ten.runRound(
        /*round_id=*/round, /*scale=*/10, /*step_size=*/10,
        /*regularization_strength=*/10);
    require(stats.disagreement_count == 1,
            "opposite-direction cycle should not agree at scale 10 alone");
    if (round > 1) {
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
      require(stats.regularization_budget == 10,
              "legacy replay should preserve its historical binary budget");
#else
      require(stats.regularization_budget ==
                  previous_regularization_budget + 10,
              "each persistent cycle round should consume another epsilon "
              "budget unit");
#endif
    }
    previous_regularization_budget = stats.regularization_budget;
  }

  auto unregularized_schedule = make_coordinator();
  long round_id = 1;
  for (int round = 0; round < 10; ++round) {
    const auto stats = unregularized_schedule.runRound(
        /*round_id=*/round_id++, /*scale=*/10, /*step_size=*/10,
        /*regularization_strength=*/0);
    require(stats.disagreement_count == 1,
            "scale 10 unregularized prefix should keep cycling");
  }

  mcpd3::PartitionWorkerRoundStats unit_stats;
  for (int round = 0; round < 100; ++round) {
    unit_stats = unregularized_schedule.runRound(
        /*round_id=*/round_id++, /*scale=*/10, /*step_size=*/1,
        /*regularization_strength=*/0);
    if (unit_stats.disagreement_count == 0) {
      break;
    }
  }
  require(unit_stats.disagreement_count == 0,
          "unit scale should resolve the cycle even without regularization");
  require(unit_stats.regularization_budget == 0,
          "forced-unregularized unit scale should report no budget");

  const auto packages = makeOppositeDirectionCyclePackages(
      /*source_terminal_capacity=*/-100, /*target_terminal_capacity=*/92);
  mcpd3::PartitionWorkerCoordinatorOptions solve_options;
  solve_options.initial_step_size = 10;
  solve_options.max_iteration_count = 20;
  solve_options.num_optimization_scales = 2;
  solve_options.patience = 99;
  solve_options.enable_group_stopping = false;
  solve_options.use_momentum = false;
  solve_options.min_step_size = 1;
  solve_options.objective_scale = 100;
  mcpd3::PartitionWorkerCoordinator solver(
      packages, makeInProcessWorkers(packages.size()), solve_options);
  const auto result = solver.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "full 10-to-1 scale schedule should resolve the cycle");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "full scale schedule should report low-scale regularized agreement");
  require(result.final_disagreement_count == 0,
          "full scale schedule should finish with agreement");
  require(result.progress_records.back().step_size == 10,
          "scaled epsilon regularization should resolve this cycle at scale 10");
  require(result.progress_records.back().disagreement_count == 0,
          "last progress record should be agreeing");
}

void lowObjectiveScaleCyclePromotesAndConverges() {
  const auto packages = makeOppositeDirectionCyclePackages(
      /*source_terminal_capacity=*/-100, /*target_terminal_capacity=*/92);
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 100;
  options.num_optimization_scales = 3;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.min_step_size = 1;
  options.objective_scale = 10;

  mcpd3::PartitionWorkerCoordinator solver(
      packages, makeInProcessWorkers(packages.size()), options);
  const auto result = solver.solve();

  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "promoted low objective scale cycle should reach agreement");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::NO_DISAGREEMENT ||
              result.stop_reason ==
                  mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "promoted low objective scale cycle should reach agreeing stop");
  require(result.final_disagreement_count == 0,
          "promoted low objective scale cycle should finish with agreement");
  require(result.objective_scale_promotion_count == 1,
          "low objective scale cycle should promote once");
  require(result.scale == 100,
          "low objective scale cycle should finish at promoted scale");
  require(result.final_regularization_budget < result.scale,
          "promoted cycle should finish under the active budget");
  require(result.progress_records.back().disagreement_count == 0,
          "last promoted cycle progress record should be agreeing");
}

void fullSolvePromotesObjectiveScaleOnOverBudget() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 2;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;
  options.max_objective_scale_promotions = 2;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{1000, 0, 10, 0, 1, 0},
                                {40, 0, 0, 0, 0, 0},
                                {400, 0, 1, 0, 1, 0}},
      std::deque<ScriptedRound>{{2000, 0, 0, 0, 0, 0},
                                {60, 1, 0, 0, 0, 0},
                                {600, 0, 0, 0, 0, 0}},
      options, &source_worker, &target_worker);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "promoted coordinator solve should reach agreement");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "promoted coordinator solve should stop on regularized agreement");
  require(result.objective_scale_promotion_count == 1,
          "coordinator should promote objective scale once");
  require(result.scale == 100,
          "coordinator result should report the promoted objective scale");
  require(result.best_lower_bound_raw == 1000,
          "agreement after promotion should certify the primary objective");
  require(result.best_certified_lower_bound_raw == 1000,
          "explicit certified lower bound should use the agreement certificate");
  require(result.best_regularized_objective_raw == 1000,
          "regularized objective diagnostic should preserve final objective");
  require(result.final_objective == 10,
          "promoted final objective should use the promoted objective scale");
  require(result.best_lower_bound > 9.999 && result.best_lower_bound < 10.001,
          "promoted raw lower bound should use the promoted objective scale");
  require(result.total_iterations == 3,
          "coordinator should retry from the promoted schedule");
  require(result.progress_records.size() == 2,
          "over-budget iteration should not be recorded as accepted progress");
  require(result.scale_results.size() == 3,
          "solve should include over-budget, promoted high-scale, and low-scale results");
  require(result.scale_results[0].status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  REGULARIZATION_BUDGET_EXCEEDED,
          "first scale should stop on over-budget regularization");
  require(source_worker->scaleFactors() == std::vector<long>{10},
          "source worker should receive objective rescale request");
  require(target_worker->scaleFactors() == std::vector<long>{10},
          "target worker should receive objective rescale request");
  require(source_worker->saturateScaleObjective() == std::vector<bool>{false},
          "strict promotion should not request worker saturation");
  require(target_worker->saturateScaleObjective() == std::vector<bool>{false},
          "strict promotion should not request worker saturation");
  require(source_worker->requests()[1].scale == 100,
          "promoted solve should restart at objective scale");
  require(source_worker->requests()[2].regularization_strength == 10,
          "promoted low-scale solve should re-enable scaled epsilon");
}

void fullSolvePromotionForwardsSaturationFlag() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 2;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;
  options.saturate_capacity_overflow = true;
  options.max_objective_scale_promotions = 1;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{1000, 0, 10, 0, 1, 0},
                                {40, 0, 0, 0, 0, 0},
                                {400, 0, 1, 0, 1, 0}},
      std::deque<ScriptedRound>{{2000, 0, 0, 0, 0, 0},
                                {60, 1, 0, 0, 0, 0},
                                {600, 0, 0, 0, 0, 0}},
      options, &source_worker, &target_worker);

  const auto result = coordinator.solve();
  require(result.objective_scale_promotion_count == 1,
          "saturated promotion test should promote once");
  require(source_worker->scaleFactors() == std::vector<long>{10},
          "source worker should receive saturated objective rescale request");
  require(target_worker->scaleFactors() == std::vector<long>{10},
          "target worker should receive saturated objective rescale request");
  require(source_worker->saturateScaleObjective() == std::vector<bool>{true},
          "source worker should receive saturation flag");
  require(target_worker->saturateScaleObjective() == std::vector<bool>{true},
          "target worker should receive saturation flag");
}

void fullSolvePromotesRepeatedlyAfterUnitScaleExhaustion() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 1;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;
  options.max_objective_scale_promotions = 2;

  std::deque<ScriptedRound> source_script;
  std::deque<ScriptedRound> target_script;
  for (int round = 0; round < 8; ++round) {
    source_script.push_back(ScriptedRound{10 + round, 0});
    target_script.push_back(
        ScriptedRound{20 + round, round == 7 ? 0 : 1});
  }
  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::move(source_script), std::move(target_script), options,
      &source_worker, &target_worker);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "repeated final-scale promotion should reach agreement");
  require(result.objective_scale_promotion_count == 2,
          "each exhausted unit scale should promote the objective");
  require(result.scale == 1000,
          "two final-scale promotions should increase the objective by 100");
  require(result.total_iterations == 8,
          "promoted schedules should retain every level through unit scale");
  const std::vector<long> expected_steps{1, 100, 10, 1,
                                         1000, 100, 10, 1};
  require(result.progress_records.size() == expected_steps.size(),
          "repeated promotion should record every extended schedule level");
  for (size_t index = 0; index < expected_steps.size(); ++index) {
    require(result.progress_records[index].step_size == expected_steps[index],
            "repeated promotion schedule step mismatch");
    require(result.progress_records[index].effective_step_size ==
                expected_steps[index],
            "promotion should use its requested coordinator restart step");
  }
  require(source_worker->scaleFactors() == std::vector<long>({10, 10}),
          "source worker should receive both objective promotions");
  require(target_worker->scaleFactors() == std::vector<long>({10, 10}),
          "target worker should receive both objective promotions");
}

void inProcessWorkerSaturatesObjectiveScaleOverflow() {
  mcpd3::PartitionPackage package;
  package.partition_id = 0;
  package.local_node_count = 1;
  package.terminal_capacities = {
      mcpd3::capacity_test_extreme_value() / 2 + 1};
  package.local_to_global = {0};

  mcpd3::InProcessPartitionWorker strict_worker;
  strict_worker.loadPartition(package);
  bool strict_threw = false;
  try {
    strict_worker.scaleObjective(2);
  } catch (const std::overflow_error &) {
    strict_threw = true;
  }
  require(strict_threw == mcpd3::capacity_is_bounded(),
          "strict objective scaling overflow should match capacity mode");

  mcpd3::InProcessPartitionWorker saturated_worker;
  saturated_worker.loadPartition(package);
  saturated_worker.scaleObjective(2, /*saturate_capacity_overflow=*/true);
  mcpd3::PartitionSolveRequest request;
  request.round_id = 1;
  request.partition_id = 0;
  (void)saturated_worker.solveRound(request);
}

void fullSolveStopsOverBudgetWhenPromotionDisabled() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;
  options.promote_objective_scale_on_overbudget = false;

  ScriptedPartitionWorker *source_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{1000, 0, 10, 0, 1, 0}},
      std::deque<ScriptedRound>{{2000, 0, 0, 0, 0, 0}}, options,
      &source_worker);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  REGULARIZATION_BUDGET_EXCEEDED,
          "disabled promotion should expose over-budget status");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZATION_BUDGET_EXCEEDED,
          "disabled promotion should expose over-budget stop reason");
  require(result.objective_scale_promotion_count == 0,
          "disabled promotion should not rescale objective");
  require(!result.has_best_lower_bound,
          "over-budget lower bound should not be accepted without promotion");
  require(result.progress_records.empty(),
          "over-budget iteration should not be recorded as accepted progress");
  require(source_worker->scaleFactors().empty(),
          "disabled promotion should not request worker rescale");
}

void inProcessCoordinatorPromotesObjectiveScaleOnOverBudget() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 20;
  options.num_optimization_scales = 3;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;
  options.max_objective_scale_promotions = 2;

  const auto packages = makeTieBreakRegularizationPackages();
  mcpd3::PartitionWorkerCoordinator coordinator(
      packages, makeInProcessWorkers(packages.size()), options);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "in-process coordinator promotion should still reach agreement");
  require(result.objective_scale_promotion_count == 1,
          "in-process coordinator should promote objective scale");
  require(result.scale == 100,
          "in-process coordinator should finish at promoted scale");
  require(coordinator.configuredInitialStepSize() == result.scale,
          "promotion must retain its restart step for a persistent solve");
  require(result.best_lower_bound_raw == 100,
          "promoted in-process solve should preserve the exact bound: got " +
              mcpd3::integer_to_string(result.best_lower_bound_raw) +
              " final_original=" +
              mcpd3::integer_to_string(result.final_objective_raw) +
              " final_certified=" +
              mcpd3::integer_to_string(
                  result.final_certified_lower_bound_raw) +
              " final_regularized=" +
              mcpd3::integer_to_string(result.final_regularized_objective_raw) +
              " budget=" +
              mcpd3::integer_to_string(result.final_regularization_budget) +
              " contribution=" + mcpd3::integer_to_string(
                  result.final_regularization_contribution));
#if !defined(MCPD_LEGACY_32BIT_DD_REPLAY)
  require(result.final_regularization_budget == 10,
          "objective promotion must leave cumulative regularization unscaled");
#endif
  require(result.final_regularization_budget < result.scale,
          "promoted in-process solve should finish under budget");
  require(result.final_disagreement_count == 0,
          "promoted in-process solve should finish with agreement");
}

void randomInitialAlphaValidationAndZeroRadiusNoop() {
  const auto packages = makeOppositeDirectionCyclePackages();

  mcpd3::PartitionWorkerCoordinatorOptions invalid_options;
  invalid_options.randomize_initial_alphas = true;
  invalid_options.initial_alpha_random_radius = -1;
  bool threw = false;
  try {
    mcpd3::PartitionWorkerCoordinator solver(
        packages, makeInProcessWorkers(packages.size()), invalid_options);
  } catch (const std::runtime_error &e) {
    threw =
        std::string(e.what()).find("initial alpha random radius") !=
        std::string::npos;
  }
  require(threw, "negative initial alpha random radius should be rejected");

  mcpd3::PartitionWorkerCoordinatorOptions zero_radius_options;
  zero_radius_options.initial_step_size = 10;
  zero_radius_options.max_iteration_count = 1;
  zero_radius_options.num_optimization_scales = 1;
  zero_radius_options.patience = 99;
  zero_radius_options.enable_group_stopping = false;
  zero_radius_options.use_momentum = false;
  zero_radius_options.min_step_size = 1;
  zero_radius_options.regularization_scheme =
      mcpd3::PartitionWorkerRegularizationScheme::NONE;
  zero_radius_options.randomize_initial_alphas = true;
  zero_radius_options.initial_alpha_random_radius = 0;
  zero_radius_options.initial_alpha_random_seed = 9;

  mcpd3::PartitionWorkerCoordinator solver(
      packages, makeInProcessWorkers(packages.size()), zero_radius_options);
  const auto result = solver.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED,
          "zero-radius random alpha init should be a no-op");
  require(result.final_disagreement_count == 1,
          "zero-radius random alpha init should preserve disagreement");
  require(result.final_regularization_budget == 0,
          "NONE regularization scheme should not use local regularization");
  require(result.progress_records.size() == 1,
          "zero-radius random alpha run should record one round");
  require(result.progress_records[0].regularization_strength == 0,
          "NONE scheme should keep local regularization disabled");
}

void randomInitialAlphaResolvesScaleTenCycle() {
  for (const int target_terminal_capacity : {2, 5, 8}) {
    const auto packages =
        makeOppositeDirectionCyclePackages(/*source_terminal_capacity=*/-10,
                                           target_terminal_capacity);

    mcpd3::PartitionWorkerCoordinatorOptions baseline_options;
    baseline_options.initial_step_size = 10;
    baseline_options.max_iteration_count = 2;
    baseline_options.num_optimization_scales = 1;
    baseline_options.patience = 99;
    baseline_options.enable_group_stopping = false;
    baseline_options.use_momentum = false;
    baseline_options.min_step_size = 1;
    baseline_options.regularization_scheme =
        mcpd3::PartitionWorkerRegularizationScheme::NONE;
    mcpd3::PartitionWorkerCoordinator baseline(
        packages, makeInProcessWorkers(packages.size()), baseline_options);
    const auto baseline_result = baseline.solve();
    require(baseline_result.status ==
                mcpd3::PartitionWorkerOptimizationStatus::
                    ITERATION_COUNT_EXCEEDED,
            "unregularized scale 10 run should keep cycling");
    require(baseline_result.final_disagreement_count == 1,
            "unregularized scale 10 baseline should remain disagreeing");
    require(baseline_result.final_regularization_budget == 0,
            "unregularized baseline should not use local regularization");
    require(baseline_result.progress_records.size() == 2,
            "baseline should record both fixed-scale rounds");
    require(baseline_result.progress_records[0].regularization_strength == 0,
            "NONE scheme should disable baseline regularization");
    require(baseline_result.progress_records[1].regularization_strength == 0,
            "NONE scheme should keep baseline regularization disabled");

    mcpd3::PartitionWorkerCoordinatorOptions missed_options =
        baseline_options;
    missed_options.randomize_initial_alphas = true;
    missed_options.initial_alpha_random_radius = 9;
    missed_options.initial_alpha_random_seed = 1;
    mcpd3::PartitionWorkerCoordinator missed_solver(
        packages, makeInProcessWorkers(packages.size()), missed_options);
    const auto missed_result = missed_solver.solve();
    require(missed_result.status ==
                mcpd3::PartitionWorkerOptimizationStatus::
                    ITERATION_COUNT_EXCEEDED,
            "a missed random initial alpha should not resolve scale 10 cycles");
    require(missed_result.final_disagreement_count == 1,
            "a missed random initial alpha should remain disagreeing");
    require(missed_result.final_regularization_budget == 0,
            "missed random initial alpha should not use local regularization");
    require(missed_result.progress_records.size() == 2,
            "missed random initial alpha should record both rounds");
    require(missed_result.progress_records[0].regularization_strength == 0,
            "missed random initial alpha should keep local regularization off");
    require(missed_result.progress_records[1].regularization_strength == 0,
            "missed random initial alpha should keep regularization disabled");

    mcpd3::PartitionWorkerCoordinatorOptions options = baseline_options;
    options.randomize_initial_alphas = true;
    options.initial_alpha_random_radius = 9;
    options.initial_alpha_random_seed = 9;

    mcpd3::PartitionWorkerCoordinator solver(
        packages, makeInProcessWorkers(packages.size()), options);
    const auto result = solver.solve();
    require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
            "seeded random initial alpha should resolve scale 10 cycles");
    require(result.stop_reason ==
                mcpd3::PartitionWorkerStopReason::NO_DISAGREEMENT,
            "random initial alpha should reach unregularized agreement");
    require(result.total_iterations == 1,
            "random initial alpha should agree on the first round");
    require(result.final_disagreement_count == 0,
            "random initial alpha should finish with agreement");
    require(result.final_regularization_budget == 0,
            "random initial alpha should not use local regularization budget");
    require(result.progress_records.size() == 1,
            "random initial alpha should record the agreeing round");
    require(result.progress_records[0].regularization_strength == 0,
            "random initial alpha test should keep local regularization off");
    require(result.progress_records[0].disagreement_count == 0,
            "random initial alpha should remove the disagreement");
  }
}

void fullSolveStopsAtIterationLimit() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 3;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}, {12, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 1}, {22, 1}}, options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED,
          "persistent disagreement should hit iteration limit");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED,
          "iteration limit stop reason mismatch");
  require(result.total_iterations == 3,
          "iteration limit should run max iterations");
  require(result.final_disagreement_count == 1,
          "final disagreement should be retained");
  require(result.best_lower_bound_raw == 34,
          "best lower bound should track improvements");
}

void fullSolveStopsOnPatienceNoProgress() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 1;
  options.enable_group_stopping = false;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}}, options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS,
          "flat lower bound should stop on patience");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::NO_LOWER_BOUND_IMPROVEMENT,
          "patience stop reason mismatch");
  require(result.total_iterations == 2,
          "patience should stop after one non-improving iteration");
  require(result.progress_records.back().iterations_since_improvement == 1,
          "progress should report iterations since improvement");
}

void fullSolveCanExhaustScaleIterationsPastPatience() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 1;
  options.enable_group_stopping = false;
  options.exhaust_scale_iterations = true;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}, {10, 0}, {10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}, {20, 1}, {20, 1}, {20, 1}},
      options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  ITERATION_COUNT_EXCEEDED,
          "exhausted flat lower bound should hit the iteration limit");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED,
          "exhausted patience test should report iteration limit");
  require(result.total_iterations == 5,
          "exhausted patience test should run the full scale budget");
}

void fullSolveRegularizedExhaustionKeepsHighScalePatience() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 1;
  options.enable_group_stopping = false;
  options.exhaust_regularized_scale_iterations = true;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}}, options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS,
          "regularized-only exhaustion should keep high-scale patience");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::NO_LOWER_BOUND_IMPROVEMENT,
          "regularized-only high-scale stop reason mismatch");
  require(result.total_iterations == 2,
          "regularized-only high-scale case should stop on patience");
}

void fullSolveRegularizedExhaustionRunsLowScaleBudget() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 1;
  options.enable_group_stopping = false;
  options.exhaust_regularized_scale_iterations = true;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}, {10, 0}, {10, 0}, {10, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}, {20, 1}, {20, 1}, {20, 1}},
      options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  ITERATION_COUNT_EXCEEDED,
          "regularized scale should exhaust to the iteration limit");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED,
          "regularized exhausted scale should report iteration limit");
  require(result.total_iterations == 5,
          "regularized exhausted scale should run the full scale budget");
}

void fullSolveStopsOnLegacyPatienceAfterDelayedImprovement() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 5;
  options.num_optimization_scales = 1;
  options.patience = 1;
  options.legacy_patience = true;
  options.enable_group_stopping = false;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {10, 0}, {11, 0}},
      std::deque<ScriptedRound>{{20, 1}, {20, 1}, {21, 1}}, options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS,
          "legacy patience should report no further progress");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::LEGACY_PATIENCE,
          "legacy patience stop reason mismatch");
  require(result.total_iterations == 3,
          "legacy patience should stop on delayed improvement");
}

void fullSolveStopsOnGroupStopping() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 25;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = true;

  std::deque<ScriptedRound> source_script;
  std::deque<ScriptedRound> target_script;
  for (int i = 0; i < 10; ++i) {
    source_script.push_back(ScriptedRound{100 + i, 0});
    target_script.push_back(ScriptedRound{0, 1});
  }
  for (int i = 0; i < 10; ++i) {
    source_script.push_back(ScriptedRound{90 + i, 0});
    target_script.push_back(ScriptedRound{0, 1});
  }
  auto coordinator = makeScriptedCoordinator(std::move(source_script),
                                             std::move(target_script), options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::NO_FURTHER_PROGRESS,
          "group stopping should report no further progress");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::GROUP_STOPPING,
          "group stopping reason mismatch");
  require(result.total_iterations == 20,
          "group stopping should evaluate two full groups");
}

void fullSolveCanExhaustScaleIterationsPastGroupStopping() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 100;
  options.max_iteration_count = 25;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = true;
  options.exhaust_scale_iterations = true;

  std::deque<ScriptedRound> source_script;
  std::deque<ScriptedRound> target_script;
  for (int i = 0; i < 10; ++i) {
    source_script.push_back(ScriptedRound{100 + i, 0});
    target_script.push_back(ScriptedRound{0, 1});
  }
  for (int i = 0; i < 15; ++i) {
    source_script.push_back(ScriptedRound{90 + i, 0});
    target_script.push_back(ScriptedRound{0, 1});
  }
  auto coordinator = makeScriptedCoordinator(std::move(source_script),
                                             std::move(target_script), options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::
                  ITERATION_COUNT_EXCEEDED,
          "exhausted group-stopping case should hit the iteration limit");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::ITERATION_COUNT_EXCEEDED,
          "exhausted group-stopping test should report iteration limit");
  require(result.total_iterations == 25,
          "exhausted group-stopping test should run the full scale budget");
}

void fullSolveRequestsRegularizationOnlyAtLowScales() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 1000;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 4;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.max_objective_scale_promotions = 0;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}, {12, 0}, {13, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 1}, {22, 1}, {23, 1}},
      options, &source_worker, &target_worker);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED,
          "persistent disagreement should exhaust all configured scales");
  require(result.progress_records.size() == 4,
          "one iteration per scale should produce four progress records");
  require(result.progress_records[0].step_size == 1000,
          "first scale step size mismatch");
  require(result.progress_records[1].step_size == 100,
          "second scale step size mismatch");
  require(result.progress_records[2].step_size == 10,
          "third scale step size mismatch");
  require(result.progress_records[3].step_size == 1,
          "fourth scale step size mismatch");
  require(source_worker->requests()[0].regularization_strength == 0,
          "scale 1000 should not be regularized");
  require(source_worker->requests()[1].regularization_strength == 0,
          "scale 100 should not be regularized");
  require(source_worker->requests()[2].regularization_strength == 10,
          "scale 10 should request scaled epsilon regularization");
  require(source_worker->requests()[3].regularization_strength == 1,
          "scale 1 should request scaled epsilon regularization");
  require(target_worker->requests()[2].regularization_strength == 10,
          "target worker should receive scale 10 regularization");
  require(target_worker->requests()[3].regularization_strength == 1,
          "target worker should receive scale 1 regularization");
}

void fullSolveNonDecimalScheduleReachesUnitStep() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 690;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 5;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.max_objective_scale_promotions = 0;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {11, 0}, {12, 0}, {13, 0}},
      std::deque<ScriptedRound>{{20, 1}, {21, 1}, {22, 1}, {23, 1}},
      options);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED,
          "persistent disagreement should exhaust the schedule");
  require(result.progress_records.size() == 4,
          "a non-decimal schedule should include one unit-step scale");
  require(result.progress_records[0].step_size == 690,
          "first non-decimal step size mismatch");
  require(result.progress_records[1].step_size == 69,
          "second non-decimal step size mismatch");
  require(result.progress_records[2].step_size == 6,
          "third non-decimal step size mismatch");
  require(result.progress_records[3].step_size == 1,
          "non-decimal schedule must finish at unit step");
}

void fullSolveContinuesAcrossScales() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 1000;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 2;
  options.patience = 99;
  options.enable_group_stopping = false;

  ScriptedPartitionWorker *source_worker = nullptr;
  ScriptedPartitionWorker *target_worker = nullptr;
  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0}, {12, 0}},
      std::deque<ScriptedRound>{{20, 1}, {22, 0}}, options, &source_worker,
      &target_worker);

  const auto result = coordinator.solve();
  require(result.status ==
              mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "second scale agreement should stop as optimal");
  require(result.scale_results.size() == 2,
          "solve should record both scale attempts");
  require(result.scale_results[0].status ==
              mcpd3::PartitionWorkerOptimizationStatus::ITERATION_COUNT_EXCEEDED,
          "first scale should hit iteration cap");
  require(result.scale_results[1].status ==
              mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "second scale should stop as optimal");
  require(result.progress_records.size() == 2,
          "two-scale solve should have two progress records");
  require(result.progress_records[0].step_size == 1000,
          "first scale step size mismatch");
  require(result.progress_records[1].step_size == 100,
          "second scale step size mismatch");
  require(source_worker->requests()[0].regularization_strength == 0,
          "first scale should not be regularized");
  require(target_worker->requests()[1].regularization_strength == 0,
          "second scale should not be regularized");
}

void coordinatorDispatchesSolveRoundsAcrossWorkersConcurrently() {
  std::atomic<int> active_solves{0};
  std::atomic<int> max_active_solves{0};
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  workers.push_back(std::make_unique<ConcurrencyProbeWorker>(
      &active_solves, &max_active_solves));
  workers.push_back(std::make_unique<ConcurrencyProbeWorker>(
      &active_solves, &max_active_solves));

  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;

  mcpd3::PartitionWorkerCoordinator coordinator(
      makeCoordinatorPackages(), std::move(workers), options);
  const auto stats =
      coordinator.runRound(/*round_id=*/1, /*scale=*/100,
                           /*step_size=*/100, /*regularization_strength=*/0);

  require(stats.disagreement_count == 0,
          "concurrency probe should return agreeing labels");
  require(max_active_solves.load() >= 2,
          "coordinator should dispatch solveRound concurrently across workers");
}

void primalDualFlowWarmStartMatchesColdPromotedSolve() {
  mcpd3::PrimalDualMinCutSolver initial(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{3, 3}, std::vector<int>{4, -2});
  initial.solve();
  const auto warm_start = initial.captureFlowWarmStart();

  mcpd3::PrimalDualMinCutSolver warmed(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{9, 9}, std::vector<int>{10, -2});
  warmed.setForceFullMinCutRecompute(true);
  warmed.restoreFlowWarmStart(warm_start);
  warmed.solve();

  mcpd3::PrimalDualMinCutSolver cold(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{9, 9}, std::vector<int>{10, -2});
  cold.solve();
  require(warmed.getMinCutValue() == cold.getMinCutValue(),
          "promoted warm solve value should match cold solve");
  require(warmed.getMinCutSolution(0) == cold.getMinCutSolution(0) &&
              warmed.getMinCutSolution(1) == cold.getMinCutSolution(1),
          "promoted warm solve labels should match cold solve");
}

void primalDualUnsolvedScaleStartsFromZeroObjective() {
  using Solver = mcpd3::PrimalDualMinCutSolver;
  alignas(Solver) unsigned char storage[sizeof(Solver)];
  std::fill(std::begin(storage), std::end(storage), 0x7f);
  Solver *solver = new (storage) Solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{1});
  try {
    solver->scaleProblem(/*scale=*/2,
                         /*saturate_capacity_overflow=*/true);
    require(solver->getMinCutValue() == 0,
            "an unsolved solver must scale a zero objective");
  } catch (...) {
    solver->~Solver();
    throw;
  }
  solver->~Solver();
}

void primalDualFlowWarmStartRejectsUnsafeReuse() {
  mcpd3::PrimalDualMinCutSolver initial(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{4, 4}, std::vector<int>{3, -2});
  initial.solve();
  const auto warm_start = initial.captureFlowWarmStart();

  requireThrows(
      [&] {
        mcpd3::PrimalDualMinCutSolver decreased(
            /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
            std::vector<int>{3, 4}, std::vector<int>{3, -2});
        decreased.restoreFlowWarmStart(warm_start);
      },
      "warm start should reject decreased arc capacities");
  requireThrows(
      [&] {
        mcpd3::PrimalDualMinCutSolver changed_topology(
            /*nnode=*/2, /*narc=*/1, std::vector<int>{1, 0},
            std::vector<int>{4, 4}, std::vector<int>{3, -2});
        changed_topology.restoreFlowWarmStart(warm_start);
      },
      "warm start should reject changed arc topology");
  requireThrows(
      [&] {
        mcpd3::PrimalDualMinCutSolver decreased_terminal(
            /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
            std::vector<int>{4, 4}, std::vector<int>{2, -2});
        decreased_terminal.restoreFlowWarmStart(warm_start);
      },
      "warm start should reject decreased terminal capacities");

  initial.setRegularizationStrength(1);
  initial.solve();
  requireThrows([&] { (void)initial.captureFlowWarmStart(); },
                "warm start should reject internally regularized state");
}

void primalDualCapacityRefreshCanResetFlowState() {
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{5, 5}, std::vector<int>{8, -8});
  solver.solve();
  solver.replaceProblemCapacities(
      std::vector<int>{7, 7}, std::vector<int>{9, -9},
      /*preserve_flow_state=*/false);
  const auto reset = solver.captureFlowWarmStart();
  require(std::all_of(reset.v_flow.begin(), reset.v_flow.end(),
                      [](const mcpd3::Capacity &flow) { return flow == 0; }),
          "capacity refresh should reset arc flow when requested");
}

void primalDualCapacityRefreshScalesFlowStateByQuantumRatio() {
  auto seeded_state = [](mcpd3::Capacity flow) {
    mcpd3::PrimalDualMinCutSolver initial(
        /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
        std::vector<int>{6, 6}, std::vector<int>{8, -8});
    initial.solve();
    auto state = initial.captureFlowWarmStart();
    state.v_flow[0] = flow;
    state.d_flow[0] = flow;
    state.d_flow[1] = -flow;
    return state;
  };

  for (const mcpd3::Capacity old_flow :
       {mcpd3::Capacity{3}, mcpd3::Capacity{-3}}) {
    mcpd3::PrimalDualMinCutSolver solver(
        /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
        std::vector<int>{6, 6}, std::vector<int>{8, -8});
    solver.restoreFlowWarmStart(seeded_state(old_flow));
    solver.replaceProblemCapacities(
        std::vector<int>{15, 15}, std::vector<int>{17, -19},
        /*preserve_flow_state=*/true,
        /*flow_scale_numerator=*/mcpd3::Objective{5},
        /*flow_scale_denominator=*/mcpd3::Objective{2});

    const auto scaled = solver.captureFlowWarmStart();
    const mcpd3::Capacity expected =
        old_flow > 0 ? mcpd3::Capacity{7} : mcpd3::Capacity{-7};
    require(scaled.v_flow == std::vector<mcpd3::Capacity>{expected},
            "capacity refresh should scale signed arc flow toward zero");
    require(scaled.d_flow ==
                std::vector<mcpd3::NodeFlow>{expected, -expected},
            "capacity refresh should recompute balance from scaled flow");
  }

  mcpd3::PrimalDualMinCutSolver reset(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{6, 6}, std::vector<int>{8, -8});
  reset.restoreFlowWarmStart(seeded_state(3));
  reset.replaceProblemCapacities(
      std::vector<int>{15, 15}, std::vector<int>{17, -19},
      /*preserve_flow_state=*/false,
      /*flow_scale_numerator=*/mcpd3::Objective{5},
      /*flow_scale_denominator=*/mcpd3::Objective{2});
  require(reset.captureFlowWarmStart().v_flow[0] == 0,
          "flow reset should take precedence over a supplied scale ratio");

  auto require_invalid_scale = [&](const mcpd3::Objective &numerator,
                                   const mcpd3::Objective &denominator) {
    mcpd3::PrimalDualMinCutSolver solver(
        /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
        std::vector<int>{6, 6}, std::vector<int>{8, -8});
    solver.restoreFlowWarmStart(seeded_state(3));
    requireThrows(
        [&] {
          solver.replaceProblemCapacities(
              std::vector<int>{15, 15}, std::vector<int>{17, -19}, true,
              numerator, denominator);
        },
        "capacity refresh should reject a non-positive flow scale ratio");
  };
  require_invalid_scale(0, 1);
  require_invalid_scale(1, 0);
  require_invalid_scale(-1, 1);
  require_invalid_scale(1, -1);

  mcpd3::PrimalDualMinCutSolver non_proportional(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{6, 6}, std::vector<int>{8, -8});
  non_proportional.restoreFlowWarmStart(seeded_state(3));
  requireThrows(
      [&] {
        non_proportional.replaceProblemCapacities(
            std::vector<int>{14, 15}, std::vector<int>{17, -19}, true,
            mcpd3::Objective{5}, mcpd3::Objective{2});
      },
      "flow scaling should reject non-proportional internal capacities");

  if (mcpd3::capacity_is_bounded()) {
    mcpd3::PrimalDualMinCutSolver overflow(
        /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
        std::vector<int>{6, 6}, std::vector<int>{8, -8});
    overflow.restoreFlowWarmStart(seeded_state(3));
    requireThrows(
        [&] {
          overflow.replaceProblemCapacities(
              std::vector<int>{15, 15}, std::vector<int>{17, -19}, true,
              std::numeric_limits<mcpd3::Objective>::max(),
              mcpd3::Objective{1});
        },
        "flow scaling should reject widened ratio arithmetic overflow");
  }
}

void incrementalCutMaintenanceMatchesFullRecomputeAcrossLabelChanges() {
  const std::vector<int> arcs{
      0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 0,
      0, 1, 1, 4, 2, 6, 3, 7};
  const std::vector<int> arc_capacities{
      3, 5, 4, 2, 6, 1, 2, 7, 5, 3, 1, 4, 7, 2, 3, 6,
      8, 1, 2, 5, 4, 3, 6, 2};
  const std::vector<int> terminal_capacities(8, 0);

  std::list<mcpd3::DualDecompositionConstraintArc> incremental_constraints;
  std::list<mcpd3::DualDecompositionConstraintArc> full_constraints;
  mcpd3::PrimalDualMinCutSolver incremental(
      /*nnode=*/8, /*narc=*/12, std::vector<int>(arcs), arc_capacities,
      terminal_capacities);
  mcpd3::PrimalDualMinCutSolver full(
      /*nnode=*/8, /*narc=*/12, std::vector<int>(arcs), arc_capacities,
      terminal_capacities);
  full.setForceFullMinCutRecompute(true);
  for (int node = 0; node < 8; ++node) {
    incremental_constraints.emplace_back(
        /*alpha=*/100, /*last_alpha=*/100, /*alpha_momentum=*/0,
        /*partition_index_source=*/0, /*partition_index_target=*/1,
        /*local_index_source=*/node, /*local_index_target=*/-1);
    full_constraints.emplace_back(
        /*alpha=*/100, /*last_alpha=*/100, /*alpha_momentum=*/0,
        /*partition_index_source=*/0, /*partition_index_target=*/1,
        /*local_index_source=*/node, /*local_index_target=*/-1);
    incremental.addSourceDualDecompositionConstraint(
        std::prev(incremental_constraints.end()));
    full.addSourceDualDecompositionConstraint(std::prev(full_constraints.end()));
  }

  std::vector<std::vector<int>> alpha_states{
      {100, 100, 100, 100, 100, 100, 100, 100},
      {-100, -100, -100, -100, -100, -100, -100, -100},
      {100, -100, 100, -100, 100, -100, 100, -100},
      {-100, 100, -100, 100, -100, 100, -100, 100},
      {-100, -100, -100, -100, 100, 100, 100, 100},
      {100, 100, 100, 100, -100, -100, -100, -100},
      {100, -100, -100, 100, 100, -100, -100, 100},
      {-100, 100, 100, -100, -100, 100, 100, -100}};
  std::mt19937 rng(20260715);
  std::uniform_int_distribution<int> alpha_dist(-150, 150);
  for (int round = 0; round < 64; ++round) {
    std::vector<int> alphas(8);
    for (int &alpha : alphas) {
      alpha = alpha_dist(rng);
    }
    alpha_states.push_back(std::move(alphas));
  }

  std::vector<int> previous_labels(8, 0);
  bool saw_adjacent_simultaneous_changes = false;
  for (int cycle = 0; cycle < 5; ++cycle) {
    for (const auto &alphas : alpha_states) {
      auto apply_alphas = [&](auto &constraints) {
        size_t node = 0;
        for (auto &constraint : constraints) {
          constraint.last_alpha = constraint.alpha;
          constraint.alpha = alphas[node++];
        }
      };
      apply_alphas(incremental_constraints);
      apply_alphas(full_constraints);

      incremental.solve();
      full.solve();
      require(incremental.getMinCutValue() == full.getMinCutValue(),
              "incremental cut value differs from full recomputation");

      std::vector<int> labels(8);
      for (int node = 0; node < 8; ++node) {
        labels[static_cast<size_t>(node)] =
            incremental.getMinCutSolution(node);
        require(labels[static_cast<size_t>(node)] ==
                    full.getMinCutSolution(node),
                "incremental cut labels differ from full recomputation");
      }
      for (size_t edge = 0; edge < arcs.size() / 2; ++edge) {
        const int source = arcs[2 * edge];
        const int target = arcs[2 * edge + 1];
        saw_adjacent_simultaneous_changes |=
            labels[static_cast<size_t>(source)] !=
                previous_labels[static_cast<size_t>(source)] &&
            labels[static_cast<size_t>(target)] !=
                previous_labels[static_cast<size_t>(target)];
      }
      previous_labels = std::move(labels);
    }
  }
  require(saw_adjacent_simultaneous_changes,
          "test did not exercise incident-edge deduplication");
}

long binaryCutValue(const std::vector<int> &labels,
                    const std::vector<int> &arcs,
                    const std::vector<int> &arc_capacities,
                    const std::vector<int> &terminal_capacities) {
  long value = 0;
  for (size_t arc = 0; arc < arcs.size() / 2; ++arc) {
    const int source = arcs[2 * arc];
    const int target = arcs[2 * arc + 1];
    if (labels[source] == 0 && labels[target] == 1) {
      value += arc_capacities[2 * arc];
    } else if (labels[source] == 1 && labels[target] == 0) {
      value += arc_capacities[2 * arc + 1];
    }
  }
  for (size_t node = 0; node < labels.size(); ++node) {
    if (labels[node] == 0 && terminal_capacities[node] < 0) {
      value -= terminal_capacities[node];
    } else if (labels[node] == 1 && terminal_capacities[node] > 0) {
      value += terminal_capacities[node];
    }
  }
  return value;
}

void canonicalCutSelectionMatchesExhaustiveLatticeExtremes() {
  std::mt19937 rng(20260712);
  std::uniform_int_distribution<int> node_count_dist(1, 5);
  std::uniform_int_distribution<int> capacity_dist(0, 4);
  std::uniform_int_distribution<int> terminal_dist(-3, 3);

  for (int trial = 0; trial < 500; ++trial) {
    const int node_count = node_count_dist(rng);
    std::vector<int> arcs;
    std::vector<int> arc_capacities;
    for (int source = 0; source < node_count; ++source) {
      for (int target = source + 1; target < node_count; ++target) {
        if ((rng() & 1U) == 0) {
          continue;
        }
        arcs.push_back(source);
        arcs.push_back(target);
        arc_capacities.push_back(capacity_dist(rng));
        arc_capacities.push_back(capacity_dist(rng));
      }
    }
    std::vector<int> terminal_capacities(static_cast<size_t>(node_count));
    for (int &capacity : terminal_capacities) {
      capacity = terminal_dist(rng);
    }

    long optimum = std::numeric_limits<long>::max();
    std::vector<int> minimum_labels(static_cast<size_t>(node_count), 1);
    std::vector<int> maximum_labels(static_cast<size_t>(node_count), 0);
    for (int mask = 0; mask < (1 << node_count); ++mask) {
      std::vector<int> labels(static_cast<size_t>(node_count));
      for (int node = 0; node < node_count; ++node) {
        labels[static_cast<size_t>(node)] = (mask >> node) & 1;
      }
      const long value = binaryCutValue(labels, arcs, arc_capacities,
                                        terminal_capacities);
      if (value < optimum) {
        optimum = value;
        minimum_labels = labels;
        maximum_labels = labels;
      } else if (value == optimum) {
        for (int node = 0; node < node_count; ++node) {
          minimum_labels[static_cast<size_t>(node)] &=
              labels[static_cast<size_t>(node)];
          maximum_labels[static_cast<size_t>(node)] |=
              labels[static_cast<size_t>(node)];
        }
      }
    }

    auto solve = [&](mcpd3::CanonicalCutSelection selection) {
      mcpd3::PrimalDualMinCutSolver solver(
          node_count, static_cast<int>(arcs.size() / 2),
          std::vector<int>(arcs), arc_capacities, terminal_capacities);
      solver.setCanonicalCutSelection(selection);
      solver.solve();
      std::vector<int> labels(static_cast<size_t>(node_count));
      for (int node = 0; node < node_count; ++node) {
        labels[static_cast<size_t>(node)] = solver.getMinCutSolution(node);
      }
      require(solver.getMinCutValue() == optimum,
              "canonical cut must preserve the primary min-cut value");
      return labels;
    };

    require(solve(mcpd3::CanonicalCutSelection::MINIMUM_LABELS) ==
                minimum_labels,
            "canonical minimum-label cut differs from exhaustive lattice");
    require(solve(mcpd3::CanonicalCutSelection::MAXIMUM_LABELS) ==
                maximum_labels,
            "canonical maximum-label cut differs from exhaustive lattice");
  }
}

void referenceGuidedCutSelectionMatchesClosestExhaustiveOptimum() {
  std::mt19937 rng(20260713);
  std::uniform_int_distribution<int> node_count_dist(1, 6);
  std::uniform_int_distribution<int> capacity_dist(0, 5);
  std::uniform_int_distribution<int> terminal_dist(-4, 4);

  for (int trial = 0; trial < 500; ++trial) {
    const int node_count = node_count_dist(rng);
    std::vector<int> arcs;
    std::vector<int> arc_capacities;
    for (int source = 0; source < node_count; ++source) {
      for (int target = source + 1; target < node_count; ++target) {
        if ((rng() & 1U) == 0) {
          continue;
        }
        arcs.push_back(source);
        arcs.push_back(target);
        arc_capacities.push_back(capacity_dist(rng));
        arc_capacities.push_back(capacity_dist(rng));
      }
    }
    std::vector<int> terminal_capacities(static_cast<size_t>(node_count));
    std::vector<int> reference(static_cast<size_t>(node_count));
    for (int node = 0; node < node_count; ++node) {
      terminal_capacities[static_cast<size_t>(node)] = terminal_dist(rng);
      reference[static_cast<size_t>(node)] = static_cast<int>(rng() & 1U);
    }

    long optimum = std::numeric_limits<long>::max();
    int closest_distance = std::numeric_limits<int>::max();
    for (int mask = 0; mask < (1 << node_count); ++mask) {
      std::vector<int> labels(static_cast<size_t>(node_count));
      int distance = 0;
      for (int node = 0; node < node_count; ++node) {
        labels[static_cast<size_t>(node)] = (mask >> node) & 1;
        distance += labels[static_cast<size_t>(node)] !=
                    reference[static_cast<size_t>(node)];
      }
      const long value = binaryCutValue(labels, arcs, arc_capacities,
                                        terminal_capacities);
      if (value < optimum) {
        optimum = value;
        closest_distance = distance;
      } else if (value == optimum) {
        closest_distance = std::min(closest_distance, distance);
      }
    }

    mcpd3::PrimalDualMinCutSolver solver(
        node_count, static_cast<int>(arcs.size() / 2),
        std::vector<int>(arcs), arc_capacities, terminal_capacities);
    solver.setReferenceCutLabels(reference);
    solver.solve();
    int actual_distance = 0;
    for (int node = 0; node < node_count; ++node) {
      actual_distance += solver.getMinCutSolution(node) !=
                         reference[static_cast<size_t>(node)];
    }
    require(solver.getMinCutValue() == optimum,
            "reference decoding changed the primary min-cut value");
    require(actual_distance == closest_distance,
            "reference decoding did not select a closest optimal cut");
  }
}

void referenceGuidedCutSelectionValidatesLabels() {
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{1, 1}, std::vector<int>{0, 0});
  requireThrows([&] { solver.setReferenceCutLabels({0}); },
                "reference decoder should reject the wrong label count");
  requireThrows([&] { solver.setReferenceCutLabels({0, 2}); },
                "reference decoder should reject non-binary labels");
  solver.setReferenceCutLabels({1, 0});
  solver.clearReferenceCutLabels();
  solver.solve();
  require(solver.getMinCutValue() == 0,
          "clearing reference labels should preserve normal solving");
}

void exactReferenceSelectionAvoidsClosureAndPreservesOptimality() {
  {
    mcpd3::PrimalDualMinCutSolver solver(
        /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
        std::vector<int>{0, 0}, std::vector<int>{0, 0});
    solver.setReferenceCutLabels({1, 0});
    solver.setReferenceCutSelection(
        mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL);
    solver.solve();
    require(solver.getMinCutSolution(0) == 1 &&
                solver.getMinCutSolution(1) == 0,
            "exact-only selection should use an exact reference cut");
    require(solver.getReferenceClosureCount() == 0,
            "exact-only selection must not run a closure decode");
  }
  {
    mcpd3::PrimalDualMinCutSolver solver(
        /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
        std::vector<int>{1});
    solver.setReferenceCutLabels({1});
    solver.setReferenceCutSelection(
        mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL);
    solver.solve();
    require(solver.getMinCutSolution(0) == 0,
            "exact-only selection must reject a non-optimal reference");
    require(solver.getMinCutValue() == 0,
            "exact-only selection changed the primary min-cut value");
    require(solver.getReferenceExactHitCount() == 0 &&
                solver.getReferenceClosureCount() == 0,
            "non-optimal exact-only reference should keep the BK cut");
  }
}

void exactReferenceSelectionMatchesExhaustiveOptimality() {
  std::mt19937 rng(20260714);
  std::uniform_int_distribution<int> node_count_dist(1, 6);
  std::uniform_int_distribution<int> capacity_dist(0, 5);
  std::uniform_int_distribution<int> terminal_dist(-4, 4);

  for (int trial = 0; trial < 500; ++trial) {
    const int node_count = node_count_dist(rng);
    std::vector<int> arcs;
    std::vector<int> arc_capacities;
    for (int source = 0; source < node_count; ++source) {
      for (int target = source + 1; target < node_count; ++target) {
        if ((rng() & 1U) == 0) {
          continue;
        }
        arcs.push_back(source);
        arcs.push_back(target);
        arc_capacities.push_back(capacity_dist(rng));
        arc_capacities.push_back(capacity_dist(rng));
      }
    }
    std::vector<int> terminal_capacities(static_cast<size_t>(node_count));
    std::vector<int> reference(static_cast<size_t>(node_count));
    for (int node = 0; node < node_count; ++node) {
      terminal_capacities[static_cast<size_t>(node)] = terminal_dist(rng);
      reference[static_cast<size_t>(node)] = static_cast<int>(rng() & 1U);
    }

    long optimum = std::numeric_limits<long>::max();
    long reference_value = 0;
    for (int mask = 0; mask < (1 << node_count); ++mask) {
      std::vector<int> labels(static_cast<size_t>(node_count));
      for (int node = 0; node < node_count; ++node) {
        labels[static_cast<size_t>(node)] = (mask >> node) & 1;
      }
      const long value = binaryCutValue(labels, arcs, arc_capacities,
                                        terminal_capacities);
      optimum = std::min(optimum, value);
      if (labels == reference) {
        reference_value = value;
      }
    }

    mcpd3::PrimalDualMinCutSolver solver(
        node_count, static_cast<int>(arcs.size() / 2),
        std::vector<int>(arcs), arc_capacities, terminal_capacities);
    solver.setReferenceCutLabels(reference);
    solver.setReferenceCutSelection(
        mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL);
    solver.solve();
    std::vector<int> actual(static_cast<size_t>(node_count));
    for (int node = 0; node < node_count; ++node) {
      actual[static_cast<size_t>(node)] = solver.getMinCutSolution(node);
    }
    require(solver.getMinCutValue() == optimum,
            "exact-only reference decoding changed the min-cut value");
    if (reference_value == optimum) {
      require(actual == reference,
              "exact-only decoding did not select an optimal reference");
    }
    require(solver.getReferenceClosureCount() == 0,
            "exact-only decoding unexpectedly ran a closure solve");
  }
}

void referenceCutCheckIntervalIsValidatedAndApplied() {
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{0});
  solver.setReferenceCutLabels({1});
  solver.setReferenceCutSelection(
      mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL);
  requireThrows([&] { solver.setReferenceCutCheckInterval(0); },
                "reference interval should reject zero");
  requireThrows([&] { solver.setReferenceCutCheckInterval(-1); },
                "reference interval should reject negative values");
  solver.setReferenceCutCheckInterval(2);
  solver.solve();
  solver.solve();
  solver.solve();
  require(solver.getReferenceDecodeCount() == 2,
          "reference interval should check the first and every second solve");
}

void exactReferenceSelectionRespectsActiveScaledEpsilon() {
  std::list<mcpd3::DualDecompositionConstraintArc> constraints;
  constraints.emplace_back(/*alpha=*/-100, /*last_alpha=*/-90,
                           /*alpha_momentum=*/0,
                           /*partition_index_source=*/0,
                           /*partition_index_target=*/1,
                           /*local_index_source=*/0,
                           /*local_index_target=*/-1);
  auto constraint = --constraints.end();
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/1, /*narc=*/0, std::vector<int>{}, std::vector<int>{},
      std::vector<int>{-100});
  solver.addSourceDualDecompositionConstraint(constraint);
  solver.setMinCutSolution(std::vector<int>{1});
  solver.setRegularizationStrength(10);
  solver.setReferenceCutLabels({1});
  solver.setReferenceCutSelection(
      mcpd3::ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL);
  solver.solve();

  require(solver.getMinCutSolution(0) == 0,
          "reference selection must preserve the active epsilon tie-break");
  require(solver.getReferenceExactHitCount() == 0,
          "a reference rejected by epsilon must not count as an exact hit");
  require(solver.getLastRegularizationBudget() == 10,
          "reference selection should preserve epsilon diagnostics");
  require(solver.getMinCutValue() == 100,
          "reference selection changed the unregularized local lower bound");
}

void dualDecompositionPropagatesReferenceCutLabels() {
  ::setenv("MCPD3_PARTITIONER", "basic", 1);
  mcpd3::DualDecompositionOptions options;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 20;
  options.initial_step_size = 1;
  options.objective_scale = 1;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.enable_group_stopping = false;
  options.track_primal_upper_bound = false;
  options.materialize_all_partition_nodes = true;
  options.reference_cut_labels = {1, 0};

  mcpd3::DualDecomposition decomposition(
      /*npartition=*/2, /*nnode=*/2, /*narc=*/1,
      std::vector<int>{0, 1}, std::vector<int>{0, 0},
      std::vector<int>{0, 0}, options);
  decomposition.solve();

  require(decomposition.getLastDisagreementCount() == 0,
          "reference-guided tied partitions should agree");
  const auto &packages = decomposition.getPartitionPackages();
  const auto snapshots = decomposition.getPartitionSnapshots();
  for (const auto &snapshot : snapshots) {
    const auto &package = packages[static_cast<size_t>(snapshot.partition_id)];
    for (size_t local = 0; local < snapshot.local_labels.size(); ++local) {
      const int global = package.local_to_global[local];
      require(snapshot.local_labels[local] ==
                  options.reference_cut_labels[static_cast<size_t>(global)],
              "local reference label does not match its global preconditioner");
    }
  }
}

void dualDecompositionPropagatesCanonicalCutSelection() {
  ::setenv("MCPD3_PARTITIONER", "basic", 1);
  mcpd3::DualDecompositionOptions options;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 20;
  options.initial_step_size = 1;
  options.objective_scale = 1;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.enable_group_stopping = false;
  options.track_primal_upper_bound = false;
  options.materialize_all_partition_nodes = true;
  options.canonical_cut_selection =
      mcpd3::CanonicalCutSelection::MAXIMUM_LABELS;

  mcpd3::DualDecomposition decomposition(
      /*npartition=*/2, /*nnode=*/2, /*narc=*/1,
      std::vector<int>{0, 1}, std::vector<int>{0, 0},
      std::vector<int>{0, 0}, options);
  decomposition.solve();

  require(decomposition.getLastDisagreementCount() == 0,
          "canonical maximum labels should agree on a tied decomposition");
  for (const auto &snapshot : decomposition.getPartitionSnapshots()) {
    for (const int label : snapshot.local_labels) {
      require(label == 1,
              "dual decomposition did not propagate canonical selection");
    }
  }
}

mcpd3::DualDecompositionOptions warmStartDdOptions() {
  mcpd3::DualDecompositionOptions options;
  options.num_optimization_scales = 1;
  options.max_iteration_count = 1000;
  options.initial_step_size = 1;
  options.min_step_size = 1;
  options.patience = 1000;
  options.exhaust_scale_iterations = true;
  options.enable_group_stopping = false;
  options.use_momentum = true;
  options.objective_scale = 1;
  options.regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::NONE;
  options.emit_partition_packages = true;
  options.construct_solvers = true;
  options.thread_count = 1;
  return options;
}

std::unique_ptr<mcpd3::DualDecomposition> makeWarmStartDecomposition(
    int capacity) {
  return std::make_unique<mcpd3::DualDecomposition>(
      /*npartition=*/2, /*nnode=*/4, /*narc=*/3,
      std::vector<int>{0, 1, 1, 2, 2, 3},
      std::vector<int>{capacity, capacity, capacity, capacity, capacity,
                       capacity},
      std::vector<int>{capacity, 0, 0, -capacity}, warmStartDdOptions());
}

void dualDecompositionWarmStartMatchesColdPromotedSolve() {
  const char *old_partitioner = std::getenv("MCPD3_PARTITIONER");
  const std::string old_value = old_partitioner ? old_partitioner : "";
  ::setenv("MCPD3_PARTITIONER", "basic", 1);

  auto initial = makeWarmStartDecomposition(3);
  initial->solve();
  const auto warm_start = initial->captureFlowWarmStart();
  require(!warm_start.partitions.empty(),
          "DD warm start should contain partition flows");
  require(!warm_start.constraints.empty(),
          "DD warm start should contain alpha state");

  auto warmed = makeWarmStartDecomposition(12);
  warmed->restoreFlowWarmStart(warm_start);
  warmed->solve();
  auto cold = makeWarmStartDecomposition(12);
  cold->solve();

  require(warmed->getLastDisagreementCount() == 0,
          "warmed promoted DD solve should agree");
  require(warmed->getLastOriginalObjectiveRaw() ==
              cold->getLastOriginalObjectiveRaw(),
          "warmed promoted DD objective should match cold solve");
  const auto warmed_partitions = warmed->getPartitionSnapshots();
  const auto cold_partitions = cold->getPartitionSnapshots();
  require(warmed_partitions.size() == cold_partitions.size(),
          "warmed partition count should match cold solve");
  for (size_t i = 0; i < warmed_partitions.size(); ++i) {
    require(warmed_partitions[i].local_labels ==
                cold_partitions[i].local_labels,
            "warmed promoted DD labels should match cold solve");
  }

  auto missing_partition = warm_start;
  missing_partition.partitions.pop_back();
  requireThrows(
      [&] {
        auto invalid = makeWarmStartDecomposition(12);
        invalid->restoreFlowWarmStart(missing_partition);
      },
      "DD warm start should reject missing partition state");

  auto wrong_constraint = warm_start;
  wrong_constraint.constraints.front().global_node_id += 1;
  requireThrows(
      [&] {
        auto invalid = makeWarmStartDecomposition(12);
        invalid->restoreFlowWarmStart(wrong_constraint);
      },
      "DD warm start should reject mismatched alpha topology");

  if (old_partitioner) {
    ::setenv("MCPD3_PARTITIONER", old_value.c_str(), 1);
  } else {
    ::unsetenv("MCPD3_PARTITIONER");
  }
}

void dualDecompositionCapacityRefreshPreservesPersistentState() {
  const char *old_partitioner = std::getenv("MCPD3_PARTITIONER");
  const std::string old_value = old_partitioner ? old_partitioner : "";
  ::setenv("MCPD3_PARTITIONER", "basic", 1);

  auto options = warmStartDdOptions();
  options.randomize_initial_alphas = true;
  options.initial_alpha_random_radius = 9;
  options.initial_alpha_random_seed = 9;
  mcpd3::DualDecomposition persistent(
      /*npartition=*/2, /*nnode=*/4, /*narc=*/3,
      std::vector<int>{0, 1, 1, 2, 2, 3},
      std::vector<int>{4, 4, 5, 5, 6, 6},
      std::vector<int>{8, 0, 0, -7}, options);
  persistent.solve();

  const auto packages_before = persistent.getPartitionPackages();
  const auto constraints_before = persistent.getConstraintSnapshots();
  require(!constraints_before.empty(),
          "capacity refresh fixture should contain a DD constraint");

  const std::vector<int> refreshed_arcs{2, 2, 9, 9, 3, 3};
  const std::vector<int> refreshed_terminals{-3, 4, 0, 11};
  persistent.replaceProblemCapacities(refreshed_arcs, refreshed_terminals);
  persistent.configureOptimizationSchedule(
      /*num_optimization_scales=*/3, /*initial_step_size=*/100,
      /*exhaust_regularized_scale_iterations=*/true);
  require(persistent.getConfiguredNumOptimizationScales() == 3,
          "schedule refresh should update the optimization scale count");
  require(persistent.getConfiguredInitialStepSize() == 100,
          "schedule refresh should update the initial step size");
  require(persistent.getConfiguredExhaustRegularizedScaleIterations(),
          "schedule refresh should update regularized exhaustion");

  const auto packages_after = persistent.getPartitionPackages();
  const auto constraints_after = persistent.getConstraintSnapshots();
  require(packages_after.size() == packages_before.size(),
          "capacity refresh should preserve partition count");
  for (size_t i = 0; i < packages_before.size(); ++i) {
    require(packages_after[i].partition_id == packages_before[i].partition_id,
            "capacity refresh should preserve partition ids");
    require(packages_after[i].arcs == packages_before[i].arcs,
            "capacity refresh should preserve local arc topology");
    require(packages_after[i].local_to_global ==
                packages_before[i].local_to_global,
            "capacity refresh should preserve local node mapping");
  }
  require(constraints_after.size() == constraints_before.size(),
          "capacity refresh should preserve constraint count");
  for (size_t i = 0; i < constraints_before.size(); ++i) {
    require(constraints_after[i].constraint_id ==
                constraints_before[i].constraint_id &&
                constraints_after[i].global_node_id ==
                    constraints_before[i].global_node_id &&
                constraints_after[i].partition_index_source ==
                    constraints_before[i].partition_index_source &&
                constraints_after[i].partition_index_target ==
                    constraints_before[i].partition_index_target &&
                constraints_after[i].local_index_source ==
                    constraints_before[i].local_index_source &&
                constraints_after[i].local_index_target ==
                    constraints_before[i].local_index_target,
            "capacity refresh should preserve constraint topology");
    require(constraints_after[i].alpha == constraints_before[i].alpha &&
                constraints_after[i].last_alpha ==
                    constraints_before[i].last_alpha &&
                constraints_after[i].alpha_momentum ==
                    constraints_before[i].alpha_momentum,
            "capacity refresh should preserve alpha and momentum state");
  }

  const auto flow_before_scale = persistent.captureFlowWarmStart();
  const auto constraints_before_scale = persistent.getConstraintSnapshots();
  std::vector<int> doubled_arcs = refreshed_arcs;
  for (int &capacity : doubled_arcs) {
    capacity *= 2;
  }
  persistent.replaceProblemCapacities(
      doubled_arcs, refreshed_terminals,
      /*preserve_alpha_state=*/true,
      /*preserve_flow_state=*/true,
      /*flow_scale_numerator=*/mcpd3::Objective{2},
      /*flow_scale_denominator=*/mcpd3::Objective{1});
  const auto flow_after_scale = persistent.captureFlowWarmStart();
  require(flow_after_scale.partitions.size() ==
              flow_before_scale.partitions.size(),
          "scaled DD flow state should preserve partition count");
  for (size_t partition = 0;
       partition < flow_before_scale.partitions.size(); ++partition) {
    const auto &before = flow_before_scale.partitions[partition];
    const auto &after = flow_after_scale.partitions[partition];
    require(after.v_flow.size() == before.v_flow.size() &&
                after.d_flow.size() == before.d_flow.size(),
            "scaled DD flow state should preserve local shapes");
    for (size_t i = 0; i < before.v_flow.size(); ++i) {
      require(after.v_flow[i] == before.v_flow[i] * 2,
              "DD refresh should forward the flow scale to every arc");
    }
    for (size_t i = 0; i < before.d_flow.size(); ++i) {
      require(after.d_flow[i] == before.d_flow[i] * 2,
              "DD refresh should rebuild every scaled local balance");
    }
  }
  const auto constraints_after_scale = persistent.getConstraintSnapshots();
  require(constraints_after_scale.size() == constraints_before_scale.size(),
          "flow scaling should preserve DD constraint count");
  for (size_t i = 0; i < constraints_before_scale.size(); ++i) {
    require(constraints_after_scale[i].alpha ==
                    constraints_before_scale[i].alpha &&
                constraints_after_scale[i].last_alpha ==
                    constraints_before_scale[i].last_alpha &&
                constraints_after_scale[i].alpha_momentum ==
                    constraints_before_scale[i].alpha_momentum,
            "flow scaling must not scale alpha or momentum state");
  }

  const auto before_rejected_scale = persistent.captureFlowWarmStart();
  std::vector<int> non_proportional_arcs = doubled_arcs;
  for (int &capacity : non_proportional_arcs) {
    capacity *= 2;
  }
  --non_proportional_arcs.back();
  requireThrows(
      [&] {
        persistent.replaceProblemCapacities(
            non_proportional_arcs, refreshed_terminals,
            /*preserve_alpha_state=*/true,
            /*preserve_flow_state=*/true,
            /*flow_scale_numerator=*/mcpd3::Objective{2},
            /*flow_scale_denominator=*/mcpd3::Objective{1});
      },
      "DD flow scaling should reject a non-proportional partition");
  const auto after_rejected_scale = persistent.captureFlowWarmStart();
  require(after_rejected_scale.partitions.size() ==
              before_rejected_scale.partitions.size(),
          "rejected DD flow scaling should preserve partition count");
  for (size_t partition = 0;
       partition < before_rejected_scale.partitions.size(); ++partition) {
    const auto &before = before_rejected_scale.partitions[partition];
    const auto &after = after_rejected_scale.partitions[partition];
    require(after.arc_capacities == before.arc_capacities &&
                after.terminal_capacities == before.terminal_capacities &&
                after.v_flow == before.v_flow && after.d_flow == before.d_flow,
            "rejected DD flow scaling must not partially mutate partitions");
  }

  persistent.solve();
  mcpd3::DualDecomposition cold(
      /*npartition=*/2, /*nnode=*/4, /*narc=*/3,
      std::vector<int>{0, 1, 1, 2, 2, 3}, doubled_arcs,
      refreshed_terminals, options);
  cold.solve();
  require(persistent.getLastDisagreementCount() == 0,
          "refreshed persistent solve should reach agreement");
  require(persistent.getLastOriginalObjectiveRaw() ==
              cold.getLastOriginalObjectiveRaw(),
          "refreshed persistent objective should match cold solve");

  persistent.replaceProblemCapacities(
      doubled_arcs, refreshed_terminals,
      /*preserve_alpha_state=*/false,
      /*preserve_flow_state=*/true);
  for (const auto &constraint : persistent.getConstraintSnapshots()) {
    require(constraint.alpha == 0 && constraint.last_alpha == 0 &&
                constraint.alpha_momentum == 0,
            "capacity refresh should reset all alpha state when requested");
  }

  requireThrows(
      [&] { persistent.replaceProblemCapacities({1, 1}, refreshed_terminals); },
      "capacity refresh should reject the wrong arc capacity count");
  requireThrows(
      [&] { persistent.replaceProblemCapacities(refreshed_arcs, {1, 2}); },
      "capacity refresh should reject the wrong terminal capacity count");
  requireThrows(
      [&] {
        persistent.replaceProblemCapacities({2, 2, -1, 9, 3, 3},
                                            refreshed_terminals);
      },
      "capacity refresh should reject negative arc capacities");

  mcpd3::DualDecomposition isolated(
      /*npartition=*/2, /*nnode=*/3, /*narc=*/1,
      std::vector<int>{0, 1}, std::vector<int>{1, 1},
      std::vector<int>{0, 0, 0}, options);
  requireThrows(
      [&] {
        isolated.replaceProblemCapacities(std::vector<int>{1, 1},
                                          std::vector<int>{0, 0, 1});
      },
      "capacity refresh should reject activating an absent isolated node");

  auto materialized_options = options;
  materialized_options.materialize_all_partition_nodes = true;
  mcpd3::DualDecomposition materialized(
      /*npartition=*/2, /*nnode=*/3, /*narc=*/1,
      std::vector<int>{0, 1}, std::vector<int>{1, 1},
      std::vector<int>{0, 0, 0}, materialized_options);
  materialized.replaceProblemCapacities(std::vector<int>{1, 1},
                                        std::vector<int>{0, 0, 1});
  materialized.solve();
  require(materialized.getLastDisagreementCount() == 0,
          "materialized isolated terminal refresh should solve");
  requireThrows(
      [&] { persistent.configureOptimizationSchedule(0, 1, false); },
      "schedule refresh should reject zero scales");
  requireThrows(
      [&] { persistent.configureOptimizationSchedule(1, 0, false); },
      "schedule refresh should reject zero initial step size");

  if (old_partitioner) {
    ::setenv("MCPD3_PARTITIONER", old_value.c_str(), 1);
  } else {
    ::unsetenv("MCPD3_PARTITIONER");
  }
}

#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
void primalDualTracksOnlyNonzeroArcFlowUpdates() {
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{5, 5}, std::vector<int>{5, -5});

  solver.solve();
  require(solver.getArcFlowUpdateCounts().empty(),
          "disabled arc-flow tracking must allocate no per-edge counters");

  mcpd3::PrimalDualMinCutSolver tracked(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{5, 5}, std::vector<int>{5, -5});
  tracked.setTrackArcFlowUpdates(true);
  tracked.solve();
  require(tracked.getArcFlowUpdateCounts() ==
              std::vector<std::uint64_t>{1},
          "the first nonzero edge-flow change must be counted once");

  tracked.solve();
  require(tracked.getArcFlowUpdateCounts() ==
              std::vector<std::uint64_t>{1},
          "an unchanged residual solve must not add a flow update");
  tracked.resetArcFlowUpdateCounts();
  require(tracked.getArcFlowUpdateCounts() ==
              std::vector<std::uint64_t>{0},
          "arc-flow heat reset must clear every edge count");
}

void dualDecompositionFlowHeatMapsBoundaryEdgesExactlyOnce() {
  const char *old_partitioner = std::getenv("MCPD3_PARTITIONER");
  const bool had_old_partitioner = old_partitioner != nullptr;
  const std::string old_value = had_old_partitioner ? old_partitioner : "";
  ::setenv("MCPD3_PARTITIONER", "basic", 1);

  mcpd3::DualDecompositionOptions options = warmStartDdOptions();
  options.track_arc_flow_updates = true;
  options.materialize_all_partition_nodes = true;
  mcpd3::DualDecomposition decomposition(
      /*npartition=*/2, /*nnode=*/4, /*narc=*/3,
      // The middle edge is deliberately reversed. It crosses the basic
      // partition boundary after mcpd3 canonicalizes its endpoints.
      std::vector<int>{0, 1, 2, 1, 2, 3},
      std::vector<int>{5, 5, 5, 5, 5, 5},
      std::vector<int>{5, 0, 0, -5}, options);
  decomposition.solve();

  const auto heat = decomposition.getArcFlowUpdateCounts();
  require(heat.size() == 3,
          "global heatmap must contain one entry per original edge");
  require(heat[1] > 0,
          "the reversed crossing edge must retain its flow activity");

  size_t crossing_edge_copies = 0;
  for (const auto &package : decomposition.getPartitionPackages()) {
    const int local_arc_count = static_cast<int>(package.arcs.size() / 2);
    for (int local_arc = 0; local_arc < local_arc_count; ++local_arc) {
      const int local_u = package.arcs[2 * local_arc];
      const int local_v = package.arcs[2 * local_arc + 1];
      const int global_u = package.local_to_global[local_u];
      const int global_v = package.local_to_global[local_v];
      if ((global_u == 1 && global_v == 2) ||
          (global_u == 2 && global_v == 1)) {
        ++crossing_edge_copies;
      }
    }
  }
  require(crossing_edge_copies == 1,
          "a boundary edge must be owned by exactly one local partition");

  decomposition.resetArcFlowUpdateCounts();
  require(decomposition.getArcFlowUpdateCounts() ==
              std::vector<std::uint64_t>({0, 0, 0}),
          "global heat reset must clear all owned local edge counters");

  if (had_old_partitioner) {
    ::setenv("MCPD3_PARTITIONER", old_value.c_str(), 1);
  } else {
    ::unsetenv("MCPD3_PARTITIONER");
  }
}
#else
void disabledBuildRejectsInspectionTelemetry() {
  mcpd3::PrimalDualMinCutSolver solver(
      /*nnode=*/2, /*narc=*/1, std::vector<int>{0, 1},
      std::vector<int>{5, 5}, std::vector<int>{5, -5});
  requireThrows([&] { solver.setTrackArcFlowUpdates(true); },
                "default build must reject per-edge flow telemetry");
  requireThrows([&] { solver.setTrackMaxflowWorkTelemetry(true); },
                "default build must reject BK work telemetry");
  require(solver.getArcFlowUpdateCounts().empty(),
          "default build must expose no per-edge counter storage");
  require(!solver.getLastSolveWorkTelemetry().enabled,
          "default build must expose disabled BK telemetry");
}
#endif

#ifdef HAVE_METIS
void metisWeightedPartitionCutsLowActivityEdges() {
  const std::vector<int> cycle{0, 1, 1, 2, 2, 3, 3, 0};
  const std::vector<std::uint64_t> weights{100, 1, 100, 1};
  const auto labels =
      mcpd3::metis_partition(/*npartition=*/2, /*narc=*/4, /*nnode=*/4,
                             cycle, &weights);
  std::uint64_t crossing_weight = 0;
  for (size_t edge = 0; edge < weights.size(); ++edge) {
    if (labels[cycle[2 * edge]] != labels[cycle[2 * edge + 1]]) {
      crossing_weight += weights[edge];
    }
  }
  require(crossing_weight == 2,
          "weighted METIS must keep high-activity cycle edges internal");

  requireThrows(
      [&] {
        const std::vector<std::uint64_t> wrong_size{1, 1};
        (void)mcpd3::metis_partition(2, 4, 4, cycle, &wrong_size);
      },
      "weighted METIS must reject a mismatched edge-weight count");
  requireThrows(
      [&] {
        const std::vector<std::uint64_t> zero_weight{100, 0, 100, 1};
        (void)mcpd3::metis_partition(2, 4, 4, cycle, &zero_weight);
      },
      "weighted METIS must reject zero edge weights");
  if (static_cast<std::uint64_t>(std::numeric_limits<idx_t>::max()) <
      std::numeric_limits<std::uint64_t>::max()) {
    requireThrows(
        [&] {
          const std::vector<std::uint64_t> oversized{
              static_cast<std::uint64_t>(
                  std::numeric_limits<idx_t>::max()) +
                  1,
              1, 1, 1};
          (void)mcpd3::metis_partition(2, 4, 4, cycle, &oversized);
        },
        "weighted METIS must reject weights outside its integer domain");
  }

  const char *old_partitioner = std::getenv("MCPD3_PARTITIONER");
  const bool had_old_partitioner = old_partitioner != nullptr;
  const std::string old_value = had_old_partitioner ? old_partitioner : "";
  ::setenv("MCPD3_PARTITIONER", "basic", 1);
  requireThrows(
      [&] {
        (void)mcpd3::configured_graph_partition(
            2, 4, 4, cycle, /*arc_capacities=*/nullptr, &weights);
      },
      "explicit edge weights must not be silently ignored by basic "
      "partitioning");
  if (had_old_partitioner) {
    ::setenv("MCPD3_PARTITIONER", old_value.c_str(), 1);
  } else {
    ::unsetenv("MCPD3_PARTITIONER");
  }
}
#endif

void nativeAndCoordinatorUseTheSameProductionDefaults() {
  const mcpd3::DualDecompositionOptions native;
  const mcpd3::PartitionWorkerCoordinatorOptions distributed;
  require(native.objective_scale == distributed.objective_scale,
          "native and coordinator neutral objective scale defaults must match");
  require(native.initial_step_size == distributed.initial_step_size,
          "native and coordinator initial step defaults must match");
  require(native.num_optimization_scales ==
              distributed.num_optimization_scales,
          "native and coordinator scale-count defaults must match");
  require(native.patience == distributed.patience,
          "native and coordinator patience defaults must match");
  require(native.use_momentum == distributed.use_momentum,
          "native and coordinator momentum defaults must match");
  require(native.enable_group_stopping == distributed.enable_group_stopping,
          "native and coordinator group stopping defaults must match");
  require(native.regularization_scheme == distributed.regularization_scheme,
          "native and coordinator regularization defaults must match");
  require(native.scaled_epsilon_max_step_size ==
              distributed.scaled_epsilon_max_step_size,
          "native and coordinator regularization cutoffs must match");
  require(native.scaled_epsilon_strength_cap ==
              distributed.scaled_epsilon_strength_cap,
          "native and coordinator regularization caps must match");
  require(native.exhaust_scale_iterations ==
              distributed.exhaust_scale_iterations &&
              native.exhaust_regularized_scale_iterations ==
                  distributed.exhaust_regularized_scale_iterations,
          "native and coordinator exhaustion defaults must match");
  require(native.max_objective_scale_promotions ==
              distributed.max_objective_scale_promotions,
          "native and coordinator promotion defaults must match");
  require(mcpd3::kDefaultObjectiveScale == 500 &&
              mcpd3::kDefaultInitialStepSize == 5000 &&
              mcpd3::kDefaultScaledEpsilonMaxStepSize == 12 &&
              mcpd3::kDefaultScaledEpsilonStrengthCap == 2,
          "the shared product policy must preserve the validated profile");
}

} // namespace

int main() {
  try {
    nativeAndCoordinatorUseTheSameProductionDefaults();
    haloDepthOneRetainsLegacyOwnerLayout();
    finiteHaloLayoutUsesBfsMembershipAndExactMultiplicity();
    haloObjectiveMultiplierUsesNodeAndEdgeLcm();
    infiniteHaloContainsEveryNodeAndEdge();
    haloLayoutRejectsInvalidInputs();
    lowerBoundCertificateSubtractsOnlyRegularizationSlack();
    solverMemoryEstimateReportsBkAndVectorBytes();
    streamingWorkerMatchesInProcessAcrossEviction();
    streamingWorkerScalesEvictedDiskPayload();
    inProcessWorkerScopesObjectiveScalingToSelectedPartitions();
    partitionWorkersUnloadSelectedPartitionsAtomically();
    inProcessPartitionWorkerMatchesDirectSolverAcrossAlphaUpdate();
    inProcessPartitionWorkerReturnsFullLabelsOnRequest();
    exportedPartitionPackagesMatchDualDecompositionRound();
    nativeProgressReportsEveryPartitionSolve();
    exportedPartitionPackagesMaterializeBoundaryDuplicates();
    disabledPartitionPackageExportPreservesNativeSolve();
    packageOnlyExportMatchesSolverBackedExport();
    haloDepthOneExportsIdenticalLegacyPackages();
    haloDepthTwoExportsScaledInducedSubproblems();
    haloDepthTwoCertifiesTheUnscaledDirectObjective();
    haloDepthTwoPromotesWhenDuplicatedNodesExhaustTheBudget();
    haloLocalObjectivesSumToScaledGlobalObjectiveExhaustively();
    haloCapacityReplacementPreservesMultiplicityScaling();
    haloPackageCoordinatorNormalizesTheObjectiveMultiplier();
#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
    haloFlowHeatAggregatesEveryLocalEdgeCopy();
#endif
    partitionWorkerCoordinatorMatchesDualDecompositionRounds();
    partitionWorkerCoordinatorMatchesDualDecompositionRegularizedRounds();
    partitionWorkerCoordinatorCapacityRefreshMatchesNativeState();
    coordinatorReconfiguresOnlyThePerSolveSchedule();
    inProcessWorkerFileBacksEveryPersistentSolverArray();
    fileBackedNativeDdMatchesResidentStateExactly();
    fileBackedNativeDdPreservesPersistentFlowExactly();
    fileBackedNativeDdRejectsInvalidStorage();
    nativeDdRecoversGlobalLabelsWithoutPartitionPackages();
    fileBackedNativeDdBacksPersistentGlobalMetadata();
    fileBackedBasicPartitioningGeneratesLabelsWithoutResidentStaging();
    packageOnlyFileBackedDdTransfersMappedPayloadWithoutCopies();
    directedStreamingDimacsMatchesGeneralReaderValue();
    dualDecompositionRegularizationSchemeControlsLowScaleStrength();
    disagreementPlateauTrackerRequiresAFullFlatWindow();
    disagreementPlateauTrackerResetsOnProgressAndScaleReset();
    disagreementPlateauOptionsValidateAndUseUnitStrength();
    disagreementPlateauModeActivatesOnARealDdPlateau();
    dualDecompositionRandomizesExportedInitialAlphas();
    dualDecompositionObjectiveScaleIsIndependentOfStepSize();
    dualDecompositionPolyakStepUsesPrimalDualGap();
    exactBoundaryCycleDetectorRejectsHashOnlyMatches();
    exactBoundaryCycleDetectorFindsRepeatedNonconstantSequence();
    exactBoundaryCycleDetectorRejectsConstantPlateaus();
    exactBoundaryCycleDetectorHonorsRequiredRepetitions();
    exactBoundaryCycleDetectorRejectsInvalidInputs();
    speculativeCycleProbeAcceptanceCoversToleranceBoundary();
    speculativeCycleReplayRetainsTheLastSolvedAlpha();
#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
    dualDecompositionReportsPerPartitionMaxflowWork();
#endif
    coherentBoundaryDirectionRequiresDominantPartitionMode();
    disagreementTransitionClassificationCoversEveryStateChange();
    dualDecompositionProbeHonorsTotalIterationBudget();
    dualDecompositionRejectsNegativeTotalIterationBudget();
    dualDecompositionUsesExplicitPartitionLabels();
    dualDecompositionRejectsInvalidExplicitPartitionLabels();
    dualDecompositionPromotesObjectiveScaleOnOverBudget();
    dualDecompositionPromotesAfterUnitScaleExhaustion();
    dualDecompositionCanRetryUnitStepWithoutMomentum();
    coordinatorScaledEpsilonControlsMatchNativeSchedule();
    coordinatorDisagreementPlateauMatchesNativePulses();
    coordinatorHonorsTotalIterationBudget();
    coordinatorCanRetryUnitStepWithoutMomentum();
    coordinatorRejectsInvalidNativeScheduleControls();
    coordinatorRunRoundReportsCertifiedRegularizedLowerBound();
    coordinatorRunRoundStrengthensCertificateAtAgreement();
    coordinatorRunRoundLeavesUnregularizedLowerBoundUnchanged();
    scaledEpsilonRegularizationUsesPreviousSinkAnchors();
    scaledEpsilonRegularizationPersistsAfterSourceTieBreak();
#if !defined(MCPD_LEGACY_32BIT_DD_REPLAY)
    scaledEpsilonRegularizationAccumulatesAcrossPersistentSinkUpdates();
#endif
    lowScaleScaledEpsilonRegularizationHandlesBoundaryTie();
    fullSolveStopsOptimalOnUnregularizedAgreement();
    fullSolveReportsProgressThroughConfiguredCallback();
    disabledProgressCallbackDoesNotFire();
    invalidProgressIntervalIsRejected();
    fullSolveReportsBestBoundUsingObjectiveScale();
    fullSolveReportsLowScaleRegularizedAgreementAsOptimal();
    fullSolveUsesScaledEpsilonRegularizationToReachAgreement();
    coordinatorRoutesMultiplePackagesToOneWorker();
    coordinatorBatchesMultiplePackagesPerWorker();
    coordinatorBalancesInitialPackagesByCpuCapacity();
    coordinatorUsesRamAsInitialAssignmentTieBreaker();
    coordinatorRejectsMalformedBatchResponses();
    inProcessWorkerBatchSolvesDistinctLoadedPartitions();
    coordinatorCollectsFinalLabelsWhenRequested();
    coordinatorSendsOnlyDirtyAlphaUpdates();
    inProcessCoordinatorSolvesWithOneWorkerOwningAllPackages();
    unitScaleResolvesOppositeDirectionCycle();
    lowObjectiveScaleCyclePromotesAndConverges();
    fullSolvePromotesObjectiveScaleOnOverBudget();
    fullSolvePromotionForwardsSaturationFlag();
    fullSolvePromotesRepeatedlyAfterUnitScaleExhaustion();
    inProcessWorkerSaturatesObjectiveScaleOverflow();
    fullSolveStopsOverBudgetWhenPromotionDisabled();
    inProcessCoordinatorPromotesObjectiveScaleOnOverBudget();
    randomInitialAlphaValidationAndZeroRadiusNoop();
    randomInitialAlphaResolvesScaleTenCycle();
    fullSolveStopsAtIterationLimit();
    fullSolveStopsOnPatienceNoProgress();
    fullSolveCanExhaustScaleIterationsPastPatience();
    fullSolveRegularizedExhaustionKeepsHighScalePatience();
    fullSolveRegularizedExhaustionRunsLowScaleBudget();
    fullSolveStopsOnLegacyPatienceAfterDelayedImprovement();
    fullSolveStopsOnGroupStopping();
    fullSolveCanExhaustScaleIterationsPastGroupStopping();
    fullSolveRequestsRegularizationOnlyAtLowScales();
    fullSolveNonDecimalScheduleReachesUnitStep();
    fullSolveContinuesAcrossScales();
    coordinatorDispatchesSolveRoundsAcrossWorkersConcurrently();
    primalDualFlowWarmStartMatchesColdPromotedSolve();
    primalDualUnsolvedScaleStartsFromZeroObjective();
    primalDualFlowWarmStartRejectsUnsafeReuse();
    primalDualCapacityRefreshCanResetFlowState();
    primalDualCapacityRefreshScalesFlowStateByQuantumRatio();
    incrementalCutMaintenanceMatchesFullRecomputeAcrossLabelChanges();
    canonicalCutSelectionMatchesExhaustiveLatticeExtremes();
    referenceGuidedCutSelectionMatchesClosestExhaustiveOptimum();
    referenceGuidedCutSelectionValidatesLabels();
    exactReferenceSelectionAvoidsClosureAndPreservesOptimality();
    exactReferenceSelectionMatchesExhaustiveOptimality();
    referenceCutCheckIntervalIsValidatedAndApplied();
    exactReferenceSelectionRespectsActiveScaledEpsilon();
    dualDecompositionPropagatesCanonicalCutSelection();
    dualDecompositionPropagatesReferenceCutLabels();
    dualDecompositionWarmStartMatchesColdPromotedSolve();
    dualDecompositionCapacityRefreshPreservesPersistentState();
#if defined(MCPD3_ENABLE_MAXFLOW_WORK_TELEMETRY)
    primalDualTracksOnlyNonzeroArcFlowUpdates();
    dualDecompositionFlowHeatMapsBoundaryEdgesExactlyOnce();
#else
    disabledBuildRejectsInspectionTelemetry();
#endif
#ifdef HAVE_METIS
    metisWeightedPartitionCutsLowActivityEdges();
#endif
  } catch (const std::exception &e) {
    std::cerr << "partition_worker_test failed: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
