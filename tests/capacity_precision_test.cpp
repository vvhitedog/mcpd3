#include <capacity.h>
#include <decomp/partition_worker.h>
#include <graph/csrgraph.h>
#include <graph/dimacs.h>
#include <maxflow/graph.h>
#include <primaldual/mcpd3.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void precisionMetadataMatchesConfiguredType() {
  require(mcpd3::capacity_mode_name() != nullptr,
          "capacity mode name must be available");
  require(mcpd3::capacity_storage_bits() == 32 ||
              mcpd3::capacity_storage_bits() == 64 ||
              mcpd3::capacity_storage_bits() == 128 ||
              mcpd3::capacity_storage_bits() == 0,
          "capacity storage bits must identify a supported mode");
}

bool tryParseCapacity(const std::string &text, mcpd3::Capacity &value) {
  return mcpd3::parse_capacity_chars(text.data(), text.data() + text.size(),
                                     value);
}

void capacityCharacterParserHandlesSignsAndBounds() {
  mcpd3::Capacity parsed = 0;
  require(tryParseCapacity("0", parsed) && parsed == 0,
          "capacity parser must accept zero");
  require(tryParseCapacity("+17", parsed) && parsed == 17,
          "capacity parser must accept an explicit positive sign");
  require(tryParseCapacity("-17", parsed) && parsed == -17,
          "capacity parser must accept negative values");

  if (mcpd3::capacity_is_bounded()) {
    const auto minimum = std::numeric_limits<mcpd3::Capacity>::min();
    const auto maximum = std::numeric_limits<mcpd3::Capacity>::max();
    require(tryParseCapacity(mcpd3::integer_to_string(minimum), parsed) &&
                parsed == minimum,
            "capacity parser must accept the configured minimum");
    require(tryParseCapacity(mcpd3::integer_to_string(maximum), parsed) &&
                parsed == maximum,
            "capacity parser must accept the configured maximum");

    const boost::multiprecision::cpp_int too_low =
        boost::multiprecision::cpp_int(mcpd3::integer_to_string(minimum)) - 1;
    const boost::multiprecision::cpp_int too_high =
        boost::multiprecision::cpp_int(mcpd3::integer_to_string(maximum)) + 1;
    require(!tryParseCapacity(too_low.convert_to<std::string>(), parsed),
            "capacity parser must reject a value below the configured range");
    require(!tryParseCapacity(too_high.convert_to<std::string>(), parsed),
            "capacity parser must reject a value above the configured range");
  }
}

void capacityCharacterParserRejectsMalformedValues() {
  mcpd3::Capacity parsed = 0;
  for (const std::string &text : {"", "+", "-", "1x", " 1", "1 "}) {
    require(!tryParseCapacity(text, parsed),
            "capacity parser must reject malformed input: '" + text + "'");
  }
}

void nativeIntegerCapacityConversionChecksRange() {
  require(mcpd3::capacity_from_integer(17L) == 17,
          "signed native integer conversion must preserve its value");
  require(mcpd3::capacity_from_integer(23UL) == 23,
          "unsigned native integer conversion must preserve its value");

#if defined(MCPD_CAPACITY_MODE_32)
  bool high_threw = false;
  bool low_threw = false;
  try {
    (void)mcpd3::capacity_from_integer(
        static_cast<long>(std::numeric_limits<mcpd3::Capacity>::max()) + 1);
  } catch (const std::overflow_error &) {
    high_threw = true;
  }
  try {
    (void)mcpd3::capacity_from_integer(
        static_cast<long>(std::numeric_limits<mcpd3::Capacity>::min()) - 1);
  } catch (const std::overflow_error &) {
    low_threw = true;
  }
  require(high_threw && low_threw,
          "native integer conversion must reject 32-bit overflow");
#elif defined(MCPD_CAPACITY_MODE_64)
  bool high_threw = false;
  try {
    (void)mcpd3::capacity_from_integer(
        std::numeric_limits<unsigned long>::max());
  } catch (const std::overflow_error &) {
    high_threw = true;
  }
  require(high_threw,
          "native integer conversion must reject unsigned 64-bit overflow");
#endif
}

void maximalCapacityRoundTripsAndSolves() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  require(capacity > 0, "test capacity must be positive");
  require(mcpd3::parse_capacity(mcpd3::integer_to_string(capacity)) == capacity,
          "capacity must round-trip through decimal text");

  using CapacityGraph =
      Graph<mcpd3::Capacity, mcpd3::Capacity, mcpd3::Objective>;
  CapacityGraph graph(2, 1);
  graph.add_node(2);
  graph.add_tweights(0, capacity, 0);
  graph.add_tweights(1, 0, capacity);
  graph.add_edge(0, 1, capacity, 0);

  const mcpd3::Objective flow = graph.maxflow();
  require(flow == mcpd3::widen_capacity(capacity),
          "maxflow must preserve the configured extreme capacity");
  require(mcpd3::parse_objective(mcpd3::integer_to_string(flow)) == flow,
          "objective must round-trip through decimal text");
}

void configuredCapacitySurvivesGraphReallocation() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  using CapacityGraph =
      Graph<mcpd3::Capacity, mcpd3::Capacity, mcpd3::Objective>;
  CapacityGraph graph(0, 0);
  graph.add_node(18);
  graph.add_tweights(0, capacity, 0);
  graph.add_tweights(17, 0, capacity);
  for (int node = 0; node < 17; ++node) {
    graph.add_edge(node, node + 1, capacity, 0);
  }
  require(graph.maxflow() == mcpd3::widen_capacity(capacity),
          "capacity must survive BK node and arc reallocation");
}

void maximalCapacitySurvivesMcpd3SolverAndWorkerStorage() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  const std::vector<int> arcs{0, 1};
  const std::vector<mcpd3::Capacity> arc_capacities{capacity, 0};
  const std::vector<mcpd3::Capacity> terminal_capacities{capacity, -capacity};

  mcpd3::PrimalDualMinCutSolver solver(
      2, 1, std::vector<int>(arcs), arc_capacities, terminal_capacities);
  solver.solve();
  require(solver.getMinCutValue() == mcpd3::widen_capacity(capacity),
          "mcpd3 solver must preserve the configured extreme capacity");

  mcpd3::PartitionPackage package;
  package.partition_id = 0;
  package.local_node_count = 2;
  package.arcs = arcs;
  package.arc_capacities = arc_capacities;
  package.terminal_capacities = terminal_capacities;
  package.local_to_global = {0, 1};

  mcpd3::InProcessPartitionWorker in_process;
  in_process.loadPartition(package);
  mcpd3::PartitionSolveRequest request;
  request.round_id = 1;
  request.partition_id = 0;
  require(in_process.solveRound(request).lower_bound ==
              mcpd3::widen_capacity(capacity),
          "in-process package must preserve the configured extreme capacity");

  mcpd3::StreamingPartitionWorker streaming;
  streaming.loadPartition(package);
  require(streaming.solveRound(request).lower_bound ==
              mcpd3::widen_capacity(capacity),
          "streaming package must preserve the configured extreme capacity");
}

void aggregateObjectiveExceedsCapacityStorage() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  const std::vector<int> arcs{0, 1, 2, 3};
  const std::vector<mcpd3::Capacity> arc_capacities{
      capacity, 0, capacity, 0};
  const std::vector<mcpd3::Capacity> terminal_capacities{
      capacity, -capacity, capacity, -capacity};

  mcpd3::PrimalDualMinCutSolver solver(
      4, 2, std::vector<int>(arcs), arc_capacities, terminal_capacities);
  solver.solve();
  const mcpd3::Objective expected = mcpd3::checked_add(
      mcpd3::widen_capacity(capacity), mcpd3::widen_capacity(capacity));
  require(solver.getMinCutValue() == expected,
          "aggregate objective must exceed capacity storage without loss");
  require(solver.getMaxFlowValue() == expected,
          "aggregate maxflow must use the widened objective domain");
}

void maximalCapacityParsesFromDimacs() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  const auto path = std::filesystem::temp_directory_path() /
                    (std::string("mcpd3_capacity_") +
                     mcpd3::capacity_mode_name() + ".max");
  {
    std::ofstream out(path, std::ios::trunc);
    require(static_cast<bool>(out), "failed to create DIMACS precision test");
    out << "p max 4 3\n"
        << "n 1 s\n"
        << "n 4 t\n"
        << "a 1 2 " << mcpd3::integer_to_string(capacity) << "\n"
        << "a 2 3 " << mcpd3::integer_to_string(capacity) << "\n"
        << "a 3 4 " << mcpd3::integer_to_string(capacity) << "\n";
  }
  auto graph = mcpd3::read_dimacs_directed_streaming(path.string());
  std::filesystem::remove(path);
  require(graph.arc_capacities ==
              std::vector<mcpd3::Capacity>({capacity, 0}),
          "DIMACS parser must preserve the configured extreme capacity");
  mcpd3::PrimalDualMinCutSolver solver(std::move(graph));
  solver.solve();
  require(solver.getMinCutValue() == mcpd3::widen_capacity(capacity),
          "DIMACS solver path must preserve the configured extreme capacity");
}

void maximalCapacitySurvivesCsrStorage() {
  const mcpd3::Capacity capacity = mcpd3::capacity_test_extreme_value();
  const auto base = std::filesystem::temp_directory_path() /
                    (std::string("mcpd3_csr_capacity_") +
                     mcpd3::capacity_mode_name());
  const auto input = base.string() + ".max";
  const auto work_dir = base.string() + "_work";
  std::filesystem::remove(input);
  std::filesystem::remove_all(work_dir);
  {
    std::ofstream out(input, std::ios::trunc);
    require(static_cast<bool>(out), "failed to create CSR precision test");
    out << "p max 4 3\n"
        << "n 1 s\n"
        << "n 4 t\n"
        << "a 1 2 " << mcpd3::integer_to_string(capacity) << "\n"
        << "a 2 3 " << mcpd3::integer_to_string(capacity) << "\n"
        << "a 3 4 " << mcpd3::integer_to_string(capacity) << "\n";
  }
  {
    auto graph = mcpd3::read_dimacs_to_csr<>(input, work_dir);
    graph.setCut({false, true});
    require(graph.getCurrentCutValue() == mcpd3::widen_capacity(capacity),
            "CSR storage must preserve the configured extreme capacity");
  }
  std::filesystem::remove(input);
  std::filesystem::remove_all(work_dir);
}

} // namespace

int main() {
  try {
    precisionMetadataMatchesConfiguredType();
    capacityCharacterParserHandlesSignsAndBounds();
    capacityCharacterParserRejectsMalformedValues();
    nativeIntegerCapacityConversionChecksRange();
    maximalCapacityRoundTripsAndSolves();
    configuredCapacitySurvivesGraphReallocation();
    maximalCapacitySurvivesMcpd3SolverAndWorkerStorage();
    aggregateObjectiveExceedsCapacityStorage();
    maximalCapacityParsesFromDimacs();
    maximalCapacitySurvivesCsrStorage();
    std::cout << "capacity_precision_test: PASS\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "capacity_precision_test: FAIL: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
