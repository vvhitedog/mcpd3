#include <graph/dimacs.h>
#include <primaldual/mcpd3.h>
#include <primaldual/pdhg_undirected.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Config {
  std::string dimacs_path;
  std::string solver = "both";
  bool print_history = false;
  mcpd3::PdhgOptions pdhg;
};

std::string requireValue(int argc, char **argv, int *index) {
  if (*index + 1 >= argc) {
    throw std::invalid_argument(std::string(argv[*index]) +
                                " requires a value");
  }
  return argv[++*index];
}

std::size_t parseSize(const std::string &name, const std::string &value) {
  if (value.empty() || value.front() == '-') {
    throw std::invalid_argument(name + " requires a nonnegative integer");
  }
  std::size_t consumed = 0;
  const unsigned long long parsed = std::stoull(value, &consumed);
  if (consumed != value.size()) {
    throw std::invalid_argument(name + " requires a nonnegative integer");
  }
  return static_cast<std::size_t>(parsed);
}

long double parseReal(const std::string &name, const std::string &value) {
  std::size_t consumed = 0;
  const long double parsed = std::stold(value, &consumed);
  if (consumed != value.size()) {
    throw std::invalid_argument(name + " requires a real number");
  }
  return parsed;
}

bool parseBoolean(const std::string &name, const std::string &value) {
  if (value == "1" || value == "true") {
    return true;
  }
  if (value == "0" || value == "false") {
    return false;
  }
  throw std::invalid_argument(name + " requires 0/1 or false/true");
}

void printUsage() {
  std::cout
      << "mcpd3_pdhg_benchmark --dimacs PATH [options]\n"
      << "  --solver existing|pdhg|both (default both)\n"
      << "  --pdhg-max-iterations N --pdhg-check-interval N\n"
      << "  --pdhg-tau X --pdhg-sigma X --pdhg-theta X\n"
      << "  --pdhg-step-size-scale X --pdhg-step-balance X\n"
      << "  --capacity-quantum X\n"
      << "  --pdhg-eps-abs X --pdhg-eps-rel X\n"
      << "  --pdhg-time-limit-seconds X\n"
      << "  --pdhg-stagnation-checks N --pdhg-stagnation-tolerance X\n"
      << "  --pdhg-use-ergodic-primal 0|1\n"
      << "  --pdhg-use-ergodic-dual 0|1\n"
      << "  --pdhg-lower-bound-safety-factor X\n"
      << "  --print-history 0|1 --verbose 0|1\n";
}

Config parseArgs(int argc, char **argv) {
  Config config;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--help") {
      printUsage();
      std::exit(EXIT_SUCCESS);
    }
    const std::string value = requireValue(argc, argv, &index);
    if (argument == "--dimacs") {
      config.dimacs_path = value;
    } else if (argument == "--solver") {
      config.solver = value == "exact" ? "existing" : value;
    } else if (argument == "--pdhg-max-iterations") {
      config.pdhg.max_iterations = parseSize(argument, value);
    } else if (argument == "--pdhg-check-interval") {
      config.pdhg.check_interval = parseSize(argument, value);
    } else if (argument == "--pdhg-tau") {
      config.pdhg.tau = parseReal(argument, value);
    } else if (argument == "--pdhg-sigma") {
      config.pdhg.sigma = parseReal(argument, value);
    } else if (argument == "--pdhg-theta") {
      config.pdhg.theta = parseReal(argument, value);
    } else if (argument == "--pdhg-step-size-scale") {
      config.pdhg.step_size_scale = parseReal(argument, value);
    } else if (argument == "--pdhg-step-balance") {
      config.pdhg.step_balance = parseReal(argument, value);
    } else if (argument == "--capacity-quantum") {
      config.pdhg.capacity_quantum = parseReal(argument, value);
    } else if (argument == "--pdhg-eps-abs") {
      config.pdhg.absolute_gap_tolerance = parseReal(argument, value);
    } else if (argument == "--pdhg-eps-rel") {
      config.pdhg.relative_gap_tolerance = parseReal(argument, value);
    } else if (argument == "--pdhg-time-limit-seconds") {
      config.pdhg.time_limit_seconds = parseReal(argument, value);
    } else if (argument == "--pdhg-stagnation-checks") {
      config.pdhg.stagnation_checks = parseSize(argument, value);
    } else if (argument == "--pdhg-stagnation-tolerance") {
      config.pdhg.stagnation_tolerance = parseReal(argument, value);
    } else if (argument == "--pdhg-use-ergodic-primal") {
      config.pdhg.use_ergodic_primal = parseBoolean(argument, value);
    } else if (argument == "--pdhg-use-ergodic-dual") {
      config.pdhg.use_ergodic_dual = parseBoolean(argument, value);
    } else if (argument == "--pdhg-lower-bound-safety-factor") {
      config.pdhg.lower_bound_safety_factor = parseReal(argument, value);
    } else if (argument == "--print-history") {
      config.print_history = parseBoolean(argument, value);
    } else if (argument == "--verbose") {
      config.pdhg.verbose = parseBoolean(argument, value);
    } else {
      throw std::invalid_argument("unknown option: " + argument);
    }
  }
  if (config.dimacs_path.empty()) {
    throw std::invalid_argument("--dimacs is required");
  }
  if (config.solver != "existing" && config.solver != "pdhg" &&
      config.solver != "both") {
    throw std::invalid_argument("--solver must be existing, pdhg, or both");
  }
  config.pdhg.record_history =
      config.print_history || config.solver == "both";
  return config;
}

long double elapsedMilliseconds(Clock::time_point start) {
  return std::chrono::duration<long double, std::milli>(Clock::now() - start)
      .count();
}

void printHistory(const mcpd3::PdhgResult &result) {
  std::cout << "pdhg_history iteration,elapsed_seconds,fractional_objective,"
               "current_cut,best_cut,current_lower,best_lower,safe_lower,"
               "gap,relative_gap,conservation_residual,x_change,p_change\n";
  for (const auto &stats : result.history) {
    std::cout << "pdhg_check " << stats.iteration << ','
              << static_cast<double>(stats.elapsed_seconds) << ','
              << static_cast<double>(stats.fractional_objective) << ','
              << mcpd3::integer_to_string(stats.current_cut_value) << ','
              << mcpd3::integer_to_string(stats.best_cut_value) << ','
              << static_cast<double>(stats.current_dual_lower_bound) << ','
              << static_cast<double>(stats.best_dual_lower_bound) << ','
              << static_cast<double>(stats.safe_lower_bound) << ','
              << static_cast<double>(stats.certified_gap) << ','
              << static_cast<double>(stats.relative_gap) << ','
              << static_cast<double>(stats.flow_conservation_residual) << ','
              << static_cast<double>(stats.x_change) << ','
              << static_cast<double>(stats.p_change) << '\n';
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    std::cout << std::unitbuf;
    const Config config = parseArgs(argc, argv);
    const auto read_start = Clock::now();
    const mcpd3::MinCutGraph graph = mcpd3::read_dimacs(config.dimacs_path);
    const long double read_ms = elapsedMilliseconds(read_start);
    std::cout << "benchmark pdhg_undirected\n"
              << "dimacs_path " << config.dimacs_path << '\n'
              << "solver " << config.solver << '\n'
              << "capacity_mode " << mcpd3::capacity_mode_name() << '\n'
              << "graph_node_count " << graph.nnode << '\n'
              << "graph_edge_count " << graph.narc << '\n'
              << "timing_read_graph_ms " << static_cast<double>(read_ms)
              << '\n';

    const bool run_existing =
        config.solver == "existing" || config.solver == "both";
    const bool run_pdhg = config.solver == "pdhg" || config.solver == "both";
    mcpd3::Objective exact_value = 0;
    long double exact_construct_ms = 0.0L;
    long double exact_solve_ms = 0.0L;
    long double exact_total_ms = 0.0L;
    std::vector<std::uint8_t> exact_source_side;
    if (run_existing) {
      const auto exact_total_start = Clock::now();
      const auto exact_construct_start = Clock::now();
      mcpd3::PrimalDualMinCutSolver exact_solver(graph);
      exact_construct_ms = elapsedMilliseconds(exact_construct_start);
      const auto exact_start = Clock::now();
      exact_solver.solve();
      exact_solve_ms = elapsedMilliseconds(exact_start);
      exact_total_ms = elapsedMilliseconds(exact_total_start);
      exact_value = exact_solver.getMinCutValue();
      exact_source_side.resize(static_cast<std::size_t>(graph.nnode));
      for (int node = 0; node < graph.nnode; ++node) {
        exact_source_side[static_cast<std::size_t>(node)] =
            exact_solver.getMinCutSolution(node) == 0 ? 1 : 0;
      }
      std::cout << "existing_cut_value "
                << mcpd3::integer_to_string(exact_value) << '\n'
                << "existing_construct_ms "
                << static_cast<double>(exact_construct_ms) << '\n'
                << "existing_solve_ms " << static_cast<double>(exact_solve_ms)
                << '\n'
                << "existing_runtime_ms " << static_cast<double>(exact_total_ms)
                << '\n';
    }

    if (run_pdhg) {
      const auto pdhg_total_start = Clock::now();
      const auto pdhg_construct_start = Clock::now();
      const mcpd3::PdhgUndirectedMinCutSolver pdhg_solver(graph);
      const long double pdhg_construct_ms =
          elapsedMilliseconds(pdhg_construct_start);
      const mcpd3::PdhgResult result = pdhg_solver.solve(config.pdhg);
      const long double pdhg_total_ms = elapsedMilliseconds(pdhg_total_start);
      std::cout
          << "pdhg_maximum_degree " << pdhg_solver.maximum_degree() << '\n'
          << "pdhg_effective_tau "
          << static_cast<double>(result.effective_tau) << '\n'
          << "pdhg_effective_sigma "
          << static_cast<double>(result.effective_sigma) << '\n'
          << "pdhg_cut_value "
          << mcpd3::integer_to_string(result.best_cut_value) << '\n'
          << "pdhg_certified_exact " << (result.certified_exact ? 1 : 0)
          << '\n'
          << "pdhg_termination_reason "
          << mcpd3::pdhg_termination_reason_name(result.termination_reason)
          << '\n'
          << "pdhg_safe_lower_bound "
          << static_cast<double>(result.safe_lower_bound) << '\n'
          << "pdhg_certified_gap "
          << static_cast<double>(result.certified_gap) << '\n'
          << "pdhg_iterations " << result.iterations << '\n'
          << "pdhg_construct_ms " << static_cast<double>(pdhg_construct_ms)
          << '\n'
          << "pdhg_solve_ms "
          << static_cast<double>(result.elapsed_seconds * 1000.0L) << '\n'
          << "pdhg_runtime_ms " << static_cast<double>(pdhg_total_ms) << '\n'
          << "pdhg_best_cut_first_iteration "
          << result.best_cut_first_iteration << '\n';

      if (run_existing) {
        const bool same_objective = result.best_cut_value == exact_value;
        const bool same_partition = result.source_side == exact_source_side;
        std::size_t true_cut_first_iteration = 0;
        long double true_cut_first_ms = 0.0L;
        bool found_true_cut = false;
        for (const auto &stats : result.history) {
          if (stats.best_cut_value == exact_value) {
            true_cut_first_iteration = stats.iteration;
            true_cut_first_ms = stats.elapsed_seconds * 1000.0L;
            found_true_cut = true;
            break;
          }
        }
        std::cout << "comparison_same_objective " << (same_objective ? 1 : 0)
                  << '\n'
                  << "comparison_same_partition "
                  << (same_partition ? 1 : 0) << '\n'
                  << "pdhg_true_cut_discovered "
                  << (found_true_cut ? 1 : 0) << '\n'
                  << "pdhg_true_cut_first_iteration "
                  << (found_true_cut
                          ? std::to_string(true_cut_first_iteration)
                          : "unavailable")
                  << '\n'
                  << "pdhg_true_cut_first_ms "
                  << (found_true_cut
                          ? std::to_string(static_cast<double>(
                                true_cut_first_ms))
                          : "unavailable")
                  << '\n';
        if (result.certified_exact && !same_objective) {
          throw std::runtime_error(
              "PDHG exact certificate disagrees with the existing solver");
        }
      }
      if (config.print_history) {
        printHistory(result);
      }
    }
  } catch (const std::exception &error) {
    std::cerr << "mcpd3_pdhg_benchmark failed: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
