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
#include <cstdlib>
#include <iostream>
#include <limits>
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
  mcpd3::DualDecompositionOptions options;
  if (argc < 2) {
    std::cout << "usage: " << argv[0]
              << " DIMACS_MAXFLOW_FILE [NUM_PARTITIONS]\n"
              << "       " << argv[0]
              << " DIMACS_MAXFLOW_FILE --partitions N "
                 "--step-policy fixed|polyak [--theta X] [--patience N] "
                 "[--max-iterations N] [--threads N]\n";
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
  auto read_graph_mem = mcpd3::get_resident_memory_usage(
      [&] { min_cut_graph_data = mcpd3::read_dimacs(dimacs_path); });
  std::cout << "read graph mem usage: " << read_graph_mem.usage_in_gb << "GB\n";

  for ( auto &cap : min_cut_graph_data.arc_capacities ) {
    cap *= 10000;
  }
  for ( auto &cap : min_cut_graph_data.terminal_capacities ) {
    cap *= 10000;
  }

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
            << " momentum=" << options.use_momentum << "\n";

  auto primal_graph2 = mcpd3::read_dimacs_to_csr(dimacs_path);

  auto primal_decoder =
      [&](const std::vector<bool> &cut, double max_lower_bound,
          const std::list<int> &disagreeing_global_indices) -> bool {
    return false;
    primal_graph2.setCut(cut);
    // auto cut_value2 = primal_graph2.getCurrentCutValue();
    // std::cout << "primal cut value2: " << cut_value2 << "\n";
    primal_graph2.narrowBandDecode(disagreeing_global_indices);
    auto cut_value2_improved = primal_graph2.getCurrentCutValue();
    std::cout << "primal cut value2 improved: " << cut_value2_improved << "\n";

    std::cout << "max_lower_bound : "
              << static_cast<long>(std::ceil(max_lower_bound)) << "\n";
    if (cut_value2_improved == static_cast<int>(std::ceil(max_lower_bound))) {
      std::cout << "breaking beause primal solution is optimal" << std::endl;
      return true;
    }
    return false;
  };

  std::unique_ptr<mcpd3::DualDecomposition> dual_decomp;
  auto dual_decomp_mem = mcpd3::get_resident_memory_usage([&] {
    dual_decomp = std::make_unique<mcpd3::DualDecomposition>(
        npartition, std::move(min_cut_graph_data), options);
  });
  std::cout << "dual decomp mem usage: " << dual_decomp_mem.usage_in_gb
            << "GB\n";
  auto microseconds =
      mcpd3::time_lambda([&] { dual_decomp->solve<true>(primal_decoder); });
  std::cout << " full loop time : " << microseconds.count() << "ms\n";

  std::cout << " total solve loop time: "
            << dual_decomp->getTotalSolveLoopTime() << "\n";
  std::cout << " best_lower_bound : " << dual_decomp->getBestLowerBound()
            << "\n";
  std::cout << " best_upper_bound : " << dual_decomp->getBestUpperBound()
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
