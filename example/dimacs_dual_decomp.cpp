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

bool parse_regularization_scheme(
    const std::string &value,
    mcpd3::DualDecompositionRegularizationScheme *scheme) {
  if (value == "scaled-epsilon" || value == "local-lexicographic") {
    *scheme =
        mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON;
    return true;
  }
  if (value == "none") {
    *scheme = mcpd3::DualDecompositionRegularizationScheme::NONE;
    return true;
  }
  return false;
}

const char *regularization_scheme_name(
    mcpd3::DualDecompositionRegularizationScheme scheme) {
  switch (scheme) {
  case mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON:
    return "scaled-epsilon";
  case mcpd3::DualDecompositionRegularizationScheme::NONE:
    return "none";
  }
  return "unknown";
}

} // namespace

int main(int argc, char *argv[]) {

  int npartition = 10;
  std::string dimacs_path;
  bool stream_symmetric_input = false;
  bool stream_directed_input = false;
  long capacity_multiplier = 1;
  mcpd3::DualDecompositionOptions options;
  if (argc < 2) {
    std::cout << "usage: " << argv[0]
              << " DIMACS_MAXFLOW_FILE [NUM_PARTITIONS]\n"
              << "       " << argv[0]
              << " DIMACS_MAXFLOW_FILE --partitions N "
                 "[--patience N] "
                 "[--max-iterations N] [--threads N] "
                 "[--regularization scaled-epsilon|none] "
                 "[--regularization-budget-limit N] "
                 "[--disable-scale-promotion] "
                 "[--max-scale-promotions N] "
                 "[--random-initial-alpha-radius N] "
                 "[--random-initial-alpha-seed N] "
                 "[--capacity-multiplier N] "
                 "[--stream-directed-input] "
                 "[--disable-primal-upper-bound] [--quiet]\n";
    std::exit(EXIT_SUCCESS);
  }
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    std::string value;
    if ((value = get_option_value(i, argc, argv, arg, "--partitions")) != "") {
      npartition = std::atoi(value.c_str());
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
                                         "--capacity-multiplier")) != "") {
      capacity_multiplier = std::atol(value.c_str());
    } else if (arg == "--disable-primal-upper-bound") {
      options.track_primal_upper_bound = false;
    } else if (arg == "--quiet") {
      options.verbose = false;
    } else if (arg == "--stream-symmetric-input") {
      stream_symmetric_input = true;
    } else if (arg == "--stream-directed-input") {
      stream_directed_input = true;
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
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--regularization")) != "") {
      if (!parse_regularization_scheme(value, &options.regularization_scheme)) {
        std::cerr << "unknown regularization scheme: " << value << "\n";
        return EXIT_FAILURE;
      }
    } else if (arg == "--disable-regularization") {
      options.regularization_scheme =
          mcpd3::DualDecompositionRegularizationScheme::NONE;
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--regularization-budget-limit")) != "") {
      options.regularization_budget_limit = std::atol(value.c_str());
    } else if (arg == "--disable-scale-promotion") {
      options.promote_objective_scale_on_overbudget = false;
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--max-scale-promotions")) != "") {
      options.max_objective_scale_promotions = std::atoi(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--random-initial-alpha-radius")) !=
               "") {
      options.randomize_initial_alphas = true;
      options.initial_alpha_random_radius = std::atol(value.c_str());
    } else if ((value = get_option_value(i, argc, argv, arg,
                                         "--random-initial-alpha-seed")) !=
               "") {
      options.initial_alpha_random_seed =
          static_cast<unsigned int>(std::strtoul(value.c_str(), nullptr, 10));
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
  if (capacity_multiplier <= 0) {
    std::cerr << "capacity multiplier must be positive\n";
    return EXIT_FAILURE;
  }
  options.objective_scale = capacity_multiplier;
  mcpd3::MinCutGraph min_cut_graph_data;
  mcpd3::ResidentMemoryUsage read_graph_mem{0, 0, 0};
  auto read_graph_microseconds = mcpd3::time_lambda([&] {
    read_graph_mem = mcpd3::get_resident_memory_usage(
        [&] {
          if (stream_symmetric_input && stream_directed_input) {
            throw std::runtime_error(
                "cannot combine symmetric and directed streaming input");
          }
          if (stream_symmetric_input) {
            min_cut_graph_data =
                mcpd3::read_dimacs_symmetric_streaming(dimacs_path);
          } else if (stream_directed_input) {
            min_cut_graph_data =
                mcpd3::read_dimacs_directed_streaming(dimacs_path);
          } else {
            min_cut_graph_data = mcpd3::read_dimacs(dimacs_path);
          }
        });
  });
  std::cout << "read graph mem usage: " << read_graph_mem.usage_in_gb << "GB\n";
  std::cout << "read_graph_time_us : " << read_graph_microseconds.count()
            << "\n";

  auto scale_graph_microseconds = mcpd3::time_lambda([&] {
    if (capacity_multiplier != 1) {
      for (auto &cap : min_cut_graph_data.arc_capacities) {
        cap *= capacity_multiplier;
      }
      for (auto &cap : min_cut_graph_data.terminal_capacities) {
        cap *= capacity_multiplier;
      }
    }
  });
  std::cout << "scale_graph_time_us : " << scale_graph_microseconds.count()
            << "\n";

  std::cout << "dd_options partitions=" << npartition
            << " step_policy=fixed"
            << " patience=" << options.patience
            << " max_iterations=" << options.max_iteration_count
            << " threads=" << options.thread_count
            << " legacy_patience=" << options.legacy_patience
            << " group_stopping=" << options.enable_group_stopping
            << " primal_upper_bound=" << options.track_primal_upper_bound
            << " verbose=" << options.verbose
            << " momentum=" << options.use_momentum
            << " stream_directed_input=" << stream_directed_input
            << " stream_symmetric_input=" << stream_symmetric_input
            << " capacity_multiplier=" << capacity_multiplier
            << " regularization="
            << regularization_scheme_name(options.regularization_scheme)
            << " regularization_budget_limit="
            << options.regularization_budget_limit
            << " promote_objective_scale_on_overbudget="
            << options.promote_objective_scale_on_overbudget
            << " max_objective_scale_promotions="
            << options.max_objective_scale_promotions
            << " randomize_initial_alphas="
            << options.randomize_initial_alphas
            << " initial_alpha_random_radius="
            << options.initial_alpha_random_radius
            << " initial_alpha_random_seed="
            << options.initial_alpha_random_seed << "\n";

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
  std::cout << " objective_scale_promotion_count : "
            << dual_decomp->getObjectiveScalePromotionCount() << "\n";
  return EXIT_SUCCESS;
}
