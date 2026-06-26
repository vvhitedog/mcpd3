// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

/**
 * A example of a dual decomposition on a DIMACS format maxflow/mincut problem.
 */

#include <io/memory.h>

#include <decomp/dualdecomp.h>
#include <graph/csrgraph.h>
#include <graph/dimacs.h>
#include <measure/timer.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

std::string get_option_value(int &i, int argc, char *argv[],
                             const std::string &arg,
                             const std::string &name) {
  const std::string prefix = name + "=";
  if (arg.rfind(prefix, 0) == 0) {
    return arg.substr(prefix.size());
  }
  if (arg == name && i + 1 < argc) {
    return argv[++i];
  }
  return "";
}

} // namespace

int main(int argc, char *argv[]) {

  int npartition = 10;
  std::string dimacs_path;
  bool stream_symmetric_input = false;
  mcpd3::DualDecompositionOptions options;
  if (argc < 2) {
    std::cout << "usage: " << argv[0]
              << " DIMACS_MAXFLOW_FILE [NUM_PARTITIONS]\n"
              << "       " << argv[0]
              << " DIMACS_MAXFLOW_FILE --partitions N "
                 "--step-policy fixed|polyak [--theta X] [--patience N] "
                 "[--max-iterations N] [--threads N] "
                 "[--disable-primal-upper-bound] [--quiet]\n";
    std::exit(EXIT_SUCCESS);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    std::string value;
    if ((value = get_option_value(i, argc, argv, arg, "--partitions")) != "") {
      npartition = std::atoi(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--step-policy")) != "") {
      if (value == "fixed") {
        options.step_policy =
            mcpd3::DualDecompositionStepPolicy::FixedScaleSchedule;
      } else if (value == "polyak") {
        options.step_policy = mcpd3::DualDecompositionStepPolicy::PolyakGap;
      } else {
        std::cerr << "unknown --step-policy: " << value << "\n";
        return EXIT_FAILURE;
      }
    } else if ((value = get_option_value(i, argc, argv, arg, "--theta")) != "") {
      options.theta = std::atof(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--patience")) != "") {
      options.patience = std::atoi(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--max-iterations")) != "") {
      options.max_iteration_count = std::atoi(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--threads")) != "") {
      options.thread_count = static_cast<size_t>(std::atoll(value.c_str()));
    } else if (arg == "--disable-primal-upper-bound") {
      options.track_primal_upper_bound = false;
    } else if (arg == "--quiet") {
      options.verbose = false;
    } else if (arg == "--stream-symmetric-input") {
      stream_symmetric_input = true;
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--min-step")) != "") {
      options.min_step_size = std::atol(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--max-step")) != "") {
      options.max_step_size = std::atol(value.c_str());
    } else if (arg == "--legacy-patience") {
      options.legacy_patience = true;
    } else if (arg == "--disable-group-stopping") {
      options.enable_group_stopping = false;
    } else if (arg == "--no-momentum") {
      options.use_momentum = false;
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "unknown option: " << arg << "\n";
      return EXIT_FAILURE;
    } else if (dimacs_path.empty()) {
      dimacs_path = arg;
    } else {
      // Backward-compatible positional partition count.
      npartition = std::atoi(arg.c_str());
    }
  }
  if (dimacs_path.empty()) {
    std::cerr << "missing DIMACS_MAXFLOW_FILE\n";
    return EXIT_FAILURE;
  }
  mcpd3::MinCutGraph min_cut_graph_data;
  mcpd3::ResidentMemoryUsage read_graph_mem{0, 0, 0};
  auto read_graph_microseconds = mcpd3::time_lambda([&] {
    read_graph_mem = mcpd3::get_resident_memory_usage(
        [&] {
          min_cut_graph_data = stream_symmetric_input
                                   ? mcpd3::read_dimacs_symmetric_streaming(
                                         dimacs_path)
                                   : mcpd3::read_dimacs(dimacs_path);
        });
  });
  std::cout << "read graph mem usage: " << read_graph_mem.usage_in_gb << "GB\n";
  std::cout << "read_graph_time_us : " << read_graph_microseconds.count()
            << "\n";

  auto scale_graph_microseconds = mcpd3::time_lambda([&] {
    for ( auto &cap : min_cut_graph_data.arc_capacities ) {
      cap *= 10000;
    }
    for ( auto &cap : min_cut_graph_data.terminal_capacities ) {
      cap *= 10000;
    }
  });
  std::cout << "scale_graph_time_us : " << scale_graph_microseconds.count()
            << "\n";

  std::cout << "dd_options partitions=" << npartition
            << " step_policy="
            << (options.step_policy ==
                        mcpd3::DualDecompositionStepPolicy::PolyakGap
                    ? "polyak"
                    : "fixed")
            << " theta=" << options.theta
            << " patience=" << options.patience
            << " max_iterations=" << options.max_iteration_count
            << " threads=" << options.thread_count
            << " legacy_patience=" << options.legacy_patience
            << " group_stopping=" << options.enable_group_stopping
            << " primal_upper_bound=" << options.track_primal_upper_bound
            << " verbose=" << options.verbose
            << " momentum=" << options.use_momentum << "\n";

  auto primal_decoder =
      [&](const std::vector<bool> &cut, double max_lower_bound,
          const std::list<int> &disagreeing_global_indices) -> bool {
    return false;
  };

  std::unique_ptr<mcpd3::DualDecomposition> dual_decomp;
  mcpd3::ResidentMemoryUsage dual_decomp_mem{0, 0, 0};
  auto dual_decomp_construct_microseconds = mcpd3::time_lambda([&] {
    dual_decomp_mem = mcpd3::get_resident_memory_usage([&] {
      dual_decomp = std::make_unique<mcpd3::DualDecomposition>(
          npartition, std::move(min_cut_graph_data), options);
    });
  });
  std::cout << "dual decomp mem usage: " << dual_decomp_mem.usage_in_gb
            << "GB\n";
  std::cout << "dual_decomp_construct_time_us : "
            << dual_decomp_construct_microseconds.count() << "\n";
  auto microseconds =
      mcpd3::time_lambda([&] { dual_decomp->solve<true>(primal_decoder); });
  std::cout << "solve_time_us : " << microseconds.count() << "\n";
  std::cout << " full loop time : " << microseconds.count() << "ms\n";

  std::cout << " total solve loop time: "
            << dual_decomp->getTotalSolveLoopTime() << "\n";
  std::cout << " iteration_count : "
            << dual_decomp->getTotalOptimizationIterations() << "\n";
  std::cout << " best_lower_bound : " << dual_decomp->getBestLowerBound()
            << "\n";
  std::cout << " best_lower_bound_raw : "
            << dual_decomp->getBestLowerBoundRaw() << "\n";
  std::cout << " best_lower_bound_unscaled : "
            << dual_decomp->getBestLowerBoundRaw() / dual_decomp->getScale()
            << "\n";
  std::cout << " best_upper_bound : " << dual_decomp->getBestUpperBound()
            << "\n";
  std::cout << " best_upper_bound_raw : "
            << dual_decomp->getBestUpperBoundRaw() << "\n";
  std::cout << " best_upper_bound_unscaled : "
            << (dual_decomp->getBestUpperBoundRaw() ==
                        std::numeric_limits<long>::max()
                    ? dual_decomp->getBestUpperBoundRaw()
                    : dual_decomp->getBestUpperBoundRaw() /
                          dual_decomp->getScale())
            << "\n";
  std::cout << " best_gap : "
            << (dual_decomp->getBestUpperBoundRaw() ==
                        std::numeric_limits<long>::max()
                    ? std::numeric_limits<double>::infinity()
                    : dual_decomp->getBestUpperBound() -
                          dual_decomp->getBestLowerBound())
            << "\n";
  std::cout << " final_disagreement_count : "
            << dual_decomp->getLastDisagreementCount() << "\n";
  std::cout << " final_disagreement_norm_sq : "
            << dual_decomp->getLastDisagreementNormSq() << "\n";
  std::cout << " final_regularization_budget_raw : "
            << dual_decomp->getLastRegularizationBudget() << "\n";
  std::cout << " final_regularization_contribution_raw : "
            << dual_decomp->getLastRegularizationContribution() << "\n";
  std::cout << " final_regularization_budget : "
            << double(dual_decomp->getLastRegularizationBudget()) /
                   dual_decomp->getScale()
            << "\n";
  std::cout << " final_regularization_contribution : "
            << double(dual_decomp->getLastRegularizationContribution()) /
                   dual_decomp->getScale()
            << "\n";
  std::cout << " final_regularization_anchor_sink_count : "
            << dual_decomp->getLastRegularizationAnchorSinkCount() << "\n";
  std::cout << " final_regularization_active_sink_count : "
            << dual_decomp->getLastRegularizationActiveSinkCount() << "\n";
  return EXIT_SUCCESS;
}
