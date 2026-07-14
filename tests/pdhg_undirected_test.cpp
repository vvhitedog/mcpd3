#include <graph/mcgraph.h>
#include <primaldual/mcpd3.h>
#include <primaldual/pdhg_undirected.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Edge = std::tuple<int, int, int>;

void require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

mcpd3::MinCutGraph makeGraph(
    int node_count, const std::vector<Edge> &edges,
    const std::vector<int> &terminal_capacities) {
  mcpd3::MinCutGraph graph;
  graph.nnode = node_count;
  graph.narc = static_cast<int>(edges.size());
  graph.terminal_capacities.reserve(terminal_capacities.size());
  for (const int capacity : terminal_capacities) {
    graph.terminal_capacities.push_back(
        mcpd3::capacity_from_integer(capacity));
  }
  for (const auto &[source, target, capacity] : edges) {
    graph.arcs.push_back(source);
    graph.arcs.push_back(target);
    graph.arc_capacities.push_back(mcpd3::capacity_from_integer(capacity));
    graph.arc_capacities.push_back(mcpd3::capacity_from_integer(capacity));
  }
  return graph;
}

mcpd3::Objective exactCutValue(const mcpd3::MinCutGraph &graph) {
  mcpd3::PrimalDualMinCutSolver solver(graph);
  solver.solve();
  return solver.getMinCutValue();
}

mcpd3::PdhgOptions testOptions() {
  mcpd3::PdhgOptions options;
  options.max_iterations = 200000;
  options.check_interval = 10;
  options.capacity_quantum = 1.0L;
  options.absolute_gap_tolerance = 1e-10L;
  options.relative_gap_tolerance = 1e-10L;
  options.use_ergodic_primal = true;
  options.use_ergodic_dual = true;
  options.record_history = true;
  return options;
}

void checkAgainstExact(const std::string &name,
                       const mcpd3::MinCutGraph &graph,
                       bool require_certificate) {
  const mcpd3::Objective exact = exactCutValue(graph);
  const long double exact_real =
      static_cast<long double>(mcpd3::integer_to_double(exact));
  const mcpd3::PdhgResult result =
      mcpd3::PdhgUndirectedMinCutSolver(graph).solve(testOptions());

  require(result.source_side.size() == static_cast<size_t>(graph.nnode),
          name + ": result partition has the wrong size");
  require(result.safe_lower_bound <= exact_real + 1e-9L,
          name + ": safe lower bound exceeds the exact value");
  require(mcpd3::widen_capacity(mcpd3::Capacity{0}) <= exact,
          name + ": exact value is unexpectedly negative");
  require(exact <= result.best_cut_value,
          name + ": best cut is below the exact value");
  require(result.certified_gap >= -1e-9L,
          name + ": certified gap is negative");
  require(!result.history.empty(), name + ": requested history is empty");
  require(result.best_cut_first_iteration <= result.iterations,
          name + ": first-best iteration is after termination");
  if (result.certified_exact) {
    require(result.best_cut_value == exact,
            name + ": exact certificate has the wrong cut value");
    require(result.termination_reason ==
                mcpd3::PdhgTerminationReason::ExactCertificate,
            name + ": exact result has the wrong termination reason");
  }
  if (require_certificate) {
    require(result.certified_exact,
            name + ": expected an exact discrete certificate");
  }
}

void requestedDeterministicGraphsRespectBounds() {
  checkAgainstExact("single edge",
                    makeGraph(2, {{0, 1, 5}}, {20, -20}), true);
  checkAgainstExact("path",
                    makeGraph(4, {{0, 1, 3}, {1, 2, 4}, {2, 3, 5}},
                              {20, 0, 0, -20}),
                    true);
  checkAgainstExact("parallel edges",
                    makeGraph(2, {{0, 1, 3}, {0, 1, 4}}, {20, -20}),
                    true);
  checkAgainstExact("disjoint paths",
                    makeGraph(4, {{0, 2, 3}, {1, 3, 5}},
                              {20, 20, -20, -20}),
                    true);
  checkAgainstExact("zero-capacity edge",
                    makeGraph(2, {{0, 1, 0}}, {20, -20}), true);
  checkAgainstExact("disconnected terminals", makeGraph(2, {}, {20, -20}),
                    true);
  checkAgainstExact("multiple optimal cuts",
                    makeGraph(2, {{0, 1, 2}}, {2, -2}), true);
}

void randomizedSmallGraphsRespectBounds() {
  std::mt19937 rng(0x50444847U);
  std::uniform_int_distribution<int> node_count_dist(2, 8);
  std::uniform_int_distribution<int> capacity_dist(0, 9);
  std::bernoulli_distribution edge_present(0.45);
  for (int trial = 0; trial < 80; ++trial) {
    const int node_count = node_count_dist(rng);
    std::vector<Edge> edges;
    for (int source = 0; source < node_count; ++source) {
      for (int target = source + 1; target < node_count; ++target) {
        if (edge_present(rng)) {
          edges.emplace_back(source, target, capacity_dist(rng));
          if ((rng() & 7U) == 0U) {
            edges.emplace_back(source, target, capacity_dist(rng));
          }
        }
      }
    }
    std::vector<int> terminals(static_cast<size_t>(node_count), 0);
    terminals.front() = 1 + capacity_dist(rng);
    terminals.back() = -(1 + capacity_dist(rng));
    if (node_count > 2 && (rng() & 1U) != 0U) {
      terminals[1] = capacity_dist(rng);
    }
    checkAgainstExact("random trial " + std::to_string(trial),
                      makeGraph(node_count, edges, terminals), false);
  }
}

void asymmetricAndMalformedGraphsAreRejected() {
  mcpd3::MinCutGraph directed = makeGraph(2, {{0, 1, 3}}, {4, -4});
  directed.arc_capacities[1] = 0;
  bool directed_threw = false;
  try {
    (void)mcpd3::PdhgUndirectedMinCutSolver(directed);
  } catch (const std::invalid_argument &) {
    directed_threw = true;
  }
  require(directed_threw, "asymmetric capacities must be rejected");

  mcpd3::MinCutGraph malformed = makeGraph(2, {{0, 1, 3}}, {4, -4});
  malformed.arcs.pop_back();
  bool malformed_threw = false;
  try {
    (void)mcpd3::PdhgUndirectedMinCutSolver(malformed);
  } catch (const std::invalid_argument &) {
    malformed_threw = true;
  }
  require(malformed_threw, "malformed storage must be rejected");
}

void terminationReasonsRemainDistinct() {
  const mcpd3::MinCutGraph graph =
      makeGraph(3, {{0, 1, 3}, {1, 2, 4}}, {20, 0, -20});

  mcpd3::PdhgOptions limited = testOptions();
  limited.max_iterations = 0;
  const auto limited_result =
      mcpd3::PdhgUndirectedMinCutSolver(graph).solve(limited);
  require(!limited_result.certified_exact,
          "iteration-limited result must not claim exactness");
  require(limited_result.termination_reason ==
              mcpd3::PdhgTerminationReason::IterationLimit,
          "iteration limit has the wrong termination reason");

  mcpd3::PdhgOptions approximate = testOptions();
  approximate.capacity_quantum = 0.0L;
  approximate.absolute_gap_tolerance = 1e9L;
  approximate.max_iterations = 10;
  approximate.check_interval = 1;
  const auto approximate_result =
      mcpd3::PdhgUndirectedMinCutSolver(graph).solve(approximate);
  require(!approximate_result.certified_exact,
          "approximate gap must not claim discrete exactness");
  require(approximate_result.termination_reason ==
              mcpd3::PdhgTerminationReason::ApproximateGap,
          "approximate gap has the wrong termination reason");

  mcpd3::PdhgOptions timed = testOptions();
  timed.capacity_quantum = 0.0L;
  timed.absolute_gap_tolerance = 0.0L;
  timed.relative_gap_tolerance = 0.0L;
  timed.time_limit_seconds = 1e-12;
  timed.max_iterations = 100;
  timed.check_interval = 1;
  const auto timed_result =
      mcpd3::PdhgUndirectedMinCutSolver(graph).solve(timed);
  require(!timed_result.certified_exact,
          "time-limited result must not claim exactness");
  require(timed_result.termination_reason ==
              mcpd3::PdhgTerminationReason::TimeLimit,
          "time limit has the wrong termination reason");

  mcpd3::PdhgOptions stagnant = testOptions();
  stagnant.capacity_quantum = 0.0L;
  stagnant.absolute_gap_tolerance = 0.0L;
  stagnant.relative_gap_tolerance = 0.0L;
  stagnant.tau = 1e-30L;
  stagnant.sigma = 1e-30L;
  stagnant.stagnation_checks = 1;
  stagnant.stagnation_tolerance = 1e-6L;
  stagnant.max_iterations = 10;
  stagnant.check_interval = 1;
  stagnant.use_ergodic_primal = false;
  stagnant.use_ergodic_dual = false;
  const auto stagnant_result =
      mcpd3::PdhgUndirectedMinCutSolver(graph).solve(stagnant);
  require(!stagnant_result.certified_exact,
          "stagnation must not claim exactness");
  require(stagnant_result.termination_reason ==
              mcpd3::PdhgTerminationReason::Stagnation,
          "stagnation has the wrong termination reason");
}

void invalidOptionsAreRejected() {
  const mcpd3::PdhgUndirectedMinCutSolver solver(
      makeGraph(2, {{0, 1, 1}}, {3, -3}));
  for (int invalid_case = 0; invalid_case < 4; ++invalid_case) {
    mcpd3::PdhgOptions options = testOptions();
    if (invalid_case == 0) {
      options.check_interval = 0;
    } else if (invalid_case == 1) {
      options.theta = 1.1L;
    } else if (invalid_case == 2) {
      options.tau = -1.0L;
    } else {
      options.capacity_quantum = -1.0L;
    }
    bool threw = false;
    try {
      (void)solver.solve(options);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    require(threw, "invalid PDHG option case was accepted");
  }
}

} // namespace

int main() {
  try {
    requestedDeterministicGraphsRespectBounds();
    randomizedSmallGraphsRespectBounds();
    asymmetricAndMalformedGraphsAreRejected();
    terminationReasonsRemainDistinct();
    invalidOptionsAreRejected();
  } catch (const std::exception &error) {
    std::cerr << "pdhg_undirected_test failed: " << error.what() << '\n';
    return 1;
  }
  std::cout << "pdhg_undirected_test passed\n";
  return 0;
}
