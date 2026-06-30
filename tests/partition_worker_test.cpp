// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <decomp/dualdecomp.h>
#include <decomp/partition_coordinator.h>
#include <decomp/partition_worker.h>
#include <graph/dimacs.h>
#include <primaldual/mcpd3.h>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

mcpd3::PartitionPackage makePackage(long alpha, long last_alpha) {
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

struct DirectSolverResult {
  long lower_bound;
  int constrained_label;
  long regularization_budget;
  long regularization_contribution;
  long regularization_anchor_sink_count;
  long regularization_active_sink_count;
};

DirectSolverResult solveDirect(mcpd3::DualDecompositionConstraintArcReference ref,
                               mcpd3::PrimalDualMinCutSolver *solver,
                               int regularization_strength) {
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
      static_cast<int>(package.arcs.size() / 2), std::vector<int>(package.arcs),
      package.arc_capacities, package.terminal_capacities);
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
  long worker_lower_bound = 0;
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

void partitionWorkerCoordinatorMatchesDualDecompositionRounds() {
  setenv("MCPD3_PARTITIONER", "basic", /*overwrite=*/1);

  auto package_source = makeTinyDualDecomposition();
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  for (size_t i = 0; i < package_source.getPartitionPackages().size(); ++i) {
    workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  }

  mcpd3::PartitionWorkerCoordinatorOptions coordinator_options;
  coordinator_options.min_step_size = 1;
  coordinator_options.max_step_size = 10000;
  coordinator_options.use_momentum = false;
  mcpd3::PartitionWorkerCoordinator coordinator(
      package_source.getPartitionPackages(), std::move(workers),
      coordinator_options);

  long best_worker_lower_bound = std::numeric_limits<long>::min();
  mcpd3::PartitionWorkerRoundStats worker_stats;
  for (long round = 1; round <= 2; ++round) {
    worker_stats = coordinator.runRound(
        /*round_id=*/round, /*scale=*/100, /*step_size=*/100,
        /*regularization_strength=*/0);
    best_worker_lower_bound =
        std::max(best_worker_lower_bound, worker_stats.lower_bound);
  }

  auto reference = makeTinyDualDecomposition();
  reference.runOptimizationScale(
      /*nstep=*/2, /*step_size=*/100, /*max_cycle_count=*/2,
      /*use_momentum=*/false);

  require(best_worker_lower_bound == reference.getBestLowerBoundRaw(),
          "coordinator best lower bound differs from DualDecomposition");
  require(worker_stats.disagreement_count ==
              reference.getLastDisagreementCount(),
          "coordinator disagreement count differs from DualDecomposition");
  require(worker_stats.disagreement_norm_sq ==
              reference.getLastDisagreementNormSq(),
          "coordinator disagreement norm differs from DualDecomposition");
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

  long randomized_alpha = 0;
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
          "over-budget regularized lower bound should not be accepted");
  require(dual_decomp.getLastRegularizationBudget() <
              dual_decomp.getScale(),
          "promoted solve should finish under the active budget limit");
  require(dual_decomp.getLastDisagreementCount() == 0,
          "promoted solve should preserve progress and reach agreement");
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
    package_ = package;
  }

  mcpd3::PartitionSolveResult solveRound(
      const mcpd3::PartitionSolveRequest &request) override {
    require(!script_.empty(), "scripted worker was called too many times");
    requests_.push_back(request);
    const auto round = script_.front();
    script_.pop_front();

    mcpd3::PartitionSolveResult result;
    result.round_id = request.round_id;
    result.partition_id = package_.partition_id;
    result.lower_bound = round.lower_bound;
    result.regularization_budget = round.regularization_budget;
    result.regularization_contribution = round.regularization_contribution;
    result.regularization_anchor_sink_count =
        round.regularization_anchor_sink_count;
    result.regularization_active_sink_count =
        round.regularization_active_sink_count;
    for (const auto &endpoint : package_.constraint_endpoints) {
      result.constrained_labels.push_back(
          mcpd3::ConstraintLabel{endpoint.constraint_id,
                                 endpoint.global_node_id,
                                 endpoint.local_index, round.label});
    }
    return result;
  }

  const std::vector<mcpd3::PartitionSolveRequest> &requests() const {
    return requests_;
  }

private:
  mcpd3::PartitionPackage package_;
  std::deque<ScriptedRound> script_;
  std::vector<mcpd3::PartitionSolveRequest> requests_;
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

std::vector<std::unique_ptr<mcpd3::PartitionWorker>> makeInProcessWorkers(
    size_t count) {
  std::vector<std::unique_ptr<mcpd3::PartitionWorker>> workers;
  for (size_t i = 0; i < count; ++i) {
    workers.push_back(std::make_unique<mcpd3::InProcessPartitionWorker>());
  }
  return workers;
}

struct OneNodeSourceSolveResult {
  long lower_bound = 0;
  int label = 0;
  long regularization_budget = 0;
  long regularization_contribution = 0;
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

void scaledEpsilonRegularizationPersistsUntilAlphaChanges() {
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
  require(solver.getLastRegularizationBudget() == 0,
          "changed alpha with previous source label should clear epsilon");
  require(solver.getMinCutSolution(/*index=*/0) == 1,
          "clearing epsilon should expose the strict sink label");
}

void lowScaleScaledEpsilonRegularizationHandlesBoundaryTie() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.min_step_size = 1;
  options.max_step_size = 10;
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
  require(record.regularization_strength == 0,
          "progress record regularization strength mismatch");
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
  require(result.progress_records.size() == 1,
          "regularized agreement should record one progress row");
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

void unitScaleResolvesOppositeDirectionCycle() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.use_momentum = false;
  options.enable_group_stopping = false;
  options.min_step_size = 1;
  options.max_step_size = 10;
  options.objective_scale = 100;

  auto make_coordinator = [&]() {
    auto packages = makeOppositeDirectionCyclePackages(
        /*source_terminal_capacity=*/-100, /*target_terminal_capacity=*/92);
    return mcpd3::PartitionWorkerCoordinator(
        packages, makeInProcessWorkers(packages.size()), options);
  };

  auto scale_ten = make_coordinator();
  for (int round = 1; round <= 4; ++round) {
    const auto stats = scale_ten.runRound(
        /*round_id=*/round, /*scale=*/10, /*step_size=*/10,
        /*regularization_strength=*/10);
    require(stats.disagreement_count == 1,
            "opposite-direction cycle should not agree at scale 10 alone");
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
  solve_options.max_step_size = 10;
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

void overBudgetScaledEpsilonAgreementRecordsDiagnostics() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 10;
  options.max_iteration_count = 2;
  options.num_optimization_scales = 1;
  options.patience = 99;
  options.enable_group_stopping = false;
  options.use_momentum = false;
  options.objective_scale = 10;

  auto coordinator = makeScriptedCoordinator(
      std::deque<ScriptedRound>{{10, 0, 10, 0, 1, 0}},
      std::deque<ScriptedRound>{{15, 0, 0, 0, 0, 0}}, options);

  const auto result = coordinator.solve();
  require(result.status == mcpd3::PartitionWorkerOptimizationStatus::OPTIMAL,
          "over-budget regularized agreement currently still stops");
  require(result.stop_reason ==
              mcpd3::PartitionWorkerStopReason::REGULARIZED_NO_DISAGREEMENT,
          "over-budget regularized agreement should report regularized stop");
  require(result.final_regularization_budget == 10,
          "over-budget diagnostics should record the strict-limit violation");
  require(result.progress_records[0].regularization_budget == 10,
          "progress should record over-budget regularization");
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
  zero_radius_options.max_step_size = 10;
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
    baseline_options.max_step_size = 10;
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

void fullSolveRequestsRegularizationOnlyAtLowScales() {
  mcpd3::PartitionWorkerCoordinatorOptions options;
  options.initial_step_size = 1000;
  options.max_iteration_count = 1;
  options.num_optimization_scales = 4;
  options.patience = 99;
  options.enable_group_stopping = false;

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

} // namespace

int main() {
  try {
    inProcessPartitionWorkerMatchesDirectSolverAcrossAlphaUpdate();
    exportedPartitionPackagesMatchDualDecompositionRound();
    partitionWorkerCoordinatorMatchesDualDecompositionRounds();
    directedStreamingDimacsMatchesGeneralReaderValue();
    dualDecompositionRegularizationSchemeControlsLowScaleStrength();
    dualDecompositionRandomizesExportedInitialAlphas();
    dualDecompositionObjectiveScaleIsIndependentOfStepSize();
    dualDecompositionPromotesObjectiveScaleOnOverBudget();
    scaledEpsilonRegularizationUsesPreviousSinkAnchors();
    scaledEpsilonRegularizationPersistsUntilAlphaChanges();
    lowScaleScaledEpsilonRegularizationHandlesBoundaryTie();
    fullSolveStopsOptimalOnUnregularizedAgreement();
    fullSolveReportsBestBoundUsingObjectiveScale();
    fullSolveReportsLowScaleRegularizedAgreementAsOptimal();
    fullSolveUsesScaledEpsilonRegularizationToReachAgreement();
    unitScaleResolvesOppositeDirectionCycle();
    overBudgetScaledEpsilonAgreementRecordsDiagnostics();
    randomInitialAlphaValidationAndZeroRadiusNoop();
    randomInitialAlphaResolvesScaleTenCycle();
    fullSolveStopsAtIterationLimit();
    fullSolveStopsOnPatienceNoProgress();
    fullSolveStopsOnLegacyPatienceAfterDelayedImprovement();
    fullSolveStopsOnGroupStopping();
    fullSolveRequestsRegularizationOnlyAtLowScales();
    fullSolveContinuesAcrossScales();
  } catch (const std::exception &e) {
    std::cerr << "partition_worker_test failed: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
