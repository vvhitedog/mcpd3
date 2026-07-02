// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#include <decomp/dualdecomp.h>
#include <graph/dimacs.h>
#include <io/memory.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

struct Timing {
  std::uint64_t total_wall_us = 0;
  std::uint64_t read_graph_wall_us = 0;
  std::uint64_t scale_graph_wall_us = 0;
  std::uint64_t construct_wall_us = 0;
  std::uint64_t solve_wall_us = 0;
};

struct MemorySnapshot {
  double vm_kb = 0;
  double rss_kb = 0;
};

struct Config {
  std::string dimacs_path;
  int partition_count = 10;
  int max_iterations = 10000;
  int schedule_levels = 5;
  long schedule_start = 10000;
  long objective_scale = 1;
  int patience = 10;
  std::size_t thread_count = 0;
  bool directed = false;
  bool symmetric_streaming = false;
  bool saturate_capacity_overflow = false;
  bool verbose = false;
  bool track_primal_upper_bound = false;
  bool emit_partition_packages = false;
  bool use_momentum = true;
  bool enable_group_stopping = true;
  bool legacy_patience = false;
  bool promote_objective_scale_on_overbudget = true;
  int max_objective_scale_promotions = 4;
  long regularization_budget_limit = 0;
  mcpd3::DualDecompositionRegularizationScheme regularization_scheme =
      mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON;
  bool randomize_initial_alphas = false;
  long initial_alpha_random_radius = 0;
  unsigned int initial_alpha_random_seed = 0;
  std::string bk_storage;
  std::string bk_mmap_dir;
  std::string bk_mmap_advise;
};

std::uint64_t elapsedUs(std::chrono::steady_clock::time_point start) {
  return static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::microseconds>(
          std::chrono::steady_clock::now() - start)
          .count());
}

MemorySnapshot memorySnapshot() {
  MemorySnapshot snapshot;
  mcpd3::process_mem_usage(snapshot.vm_kb, snapshot.rss_kb);
  return snapshot;
}

long parseLong(const std::string &value, const std::string &name) {
  char *end = nullptr;
  const long parsed = std::strtol(value.c_str(), &end, 10);
  if (end == value.c_str() || *end != '\0') {
    throw std::runtime_error("invalid integer for " + name + ": " + value);
  }
  return parsed;
}

int parseInt(const std::string &value, const std::string &name) {
  const long parsed = parseLong(value, name);
  if (parsed < std::numeric_limits<int>::min() ||
      parsed > std::numeric_limits<int>::max()) {
    throw std::runtime_error("integer out of range for " + name + ": " +
                             value);
  }
  return static_cast<int>(parsed);
}

std::size_t parseSize(const std::string &value, const std::string &name) {
  const long parsed = parseLong(value, name);
  if (parsed < 0) {
    throw std::runtime_error(name + " must be nonnegative");
  }
  return static_cast<std::size_t>(parsed);
}

bool parseRegularizationScheme(
    const std::string &value,
    mcpd3::DualDecompositionRegularizationScheme *scheme) {
  if (value == "scaled-epsilon" || value == "local-lexicographic") {
    *scheme = mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON;
    return true;
  }
  if (value == "none") {
    *scheme = mcpd3::DualDecompositionRegularizationScheme::NONE;
    return true;
  }
  return false;
}

const char *regularizationSchemeName(
    mcpd3::DualDecompositionRegularizationScheme scheme) {
  switch (scheme) {
  case mcpd3::DualDecompositionRegularizationScheme::SCALED_EPSILON:
    return "scaled-epsilon";
  case mcpd3::DualDecompositionRegularizationScheme::NONE:
    return "none";
  }
  return "unknown";
}

void printUsage(const char *program) {
  std::cerr
      << "usage: " << program
      << " DIMACS [--directed] [--symmetric-streaming]\n"
      << "       [--partitions N] [--max-iterations N]\n"
      << "       [--schedule-start N] [--schedule-levels N]\n"
      << "       [--objective-scale N] [--threads N] [--patience N]\n"
      << "       [--regularization scaled-epsilon|none]\n"
      << "       [--regularization-budget-limit N]\n"
      << "       [--disable-scale-promotion] [--max-scale-promotions N]\n"
      << "       [--random-initial-alpha-radius N]\n"
      << "       [--random-initial-alpha-seed N]\n"
      << "       [--track-primal-upper-bound] [--emit-partition-packages]\n"
      << "       [--saturate-capacity-overflow] [--verbose]\n"
      << "       [--bk-storage malloc|file_mmap|anon_mmap]\n"
      << "       [--bk-mmap-dir DIR] [--bk-mmap-advise ADVISE]\n";
}

Config parseArgs(int argc, char **argv) {
  Config config;
  if (argc < 2) {
    printUsage(argv[0]);
    throw std::runtime_error("missing DIMACS path");
  }
  config.dimacs_path = argv[1];
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto requireValue = [&](const std::string &name) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error(name + " requires a value");
      }
      return argv[++i];
    };
    if (arg == "--partitions") {
      config.partition_count = parseInt(requireValue(arg), arg);
    } else if (arg == "--max-iterations") {
      config.max_iterations = parseInt(requireValue(arg), arg);
    } else if (arg == "--schedule-levels" || arg == "--num-scales") {
      config.schedule_levels = parseInt(requireValue(arg), arg);
    } else if (arg == "--schedule-start" || arg == "--initial-step") {
      config.schedule_start = parseLong(requireValue(arg), arg);
    } else if (arg == "--objective-scale" ||
               arg == "--capacity-multiplier") {
      config.objective_scale = parseLong(requireValue(arg), arg);
    } else if (arg == "--threads") {
      config.thread_count = parseSize(requireValue(arg), arg);
    } else if (arg == "--patience") {
      config.patience = parseInt(requireValue(arg), arg);
    } else if (arg == "--directed" || arg == "--stream-directed-input") {
      config.directed = true;
    } else if (arg == "--symmetric-streaming" ||
               arg == "--stream-symmetric-input") {
      config.symmetric_streaming = true;
    } else if (arg == "--saturate-capacity-overflow") {
      config.saturate_capacity_overflow = true;
    } else if (arg == "--verbose") {
      config.verbose = true;
    } else if (arg == "--quiet") {
      config.verbose = false;
    } else if (arg == "--track-primal-upper-bound") {
      config.track_primal_upper_bound = true;
    } else if (arg == "--disable-primal-upper-bound") {
      config.track_primal_upper_bound = false;
    } else if (arg == "--emit-partition-packages") {
      config.emit_partition_packages = true;
    } else if (arg == "--no-momentum") {
      config.use_momentum = false;
    } else if (arg == "--disable-group-stopping") {
      config.enable_group_stopping = false;
    } else if (arg == "--legacy-patience") {
      config.legacy_patience = true;
    } else if (arg == "--regularization") {
      if (!parseRegularizationScheme(requireValue(arg),
                                     &config.regularization_scheme)) {
        throw std::runtime_error("unknown regularization scheme");
      }
    } else if (arg == "--disable-regularization") {
      config.regularization_scheme =
          mcpd3::DualDecompositionRegularizationScheme::NONE;
    } else if (arg == "--regularization-budget-limit") {
      config.regularization_budget_limit = parseLong(requireValue(arg), arg);
    } else if (arg == "--disable-scale-promotion") {
      config.promote_objective_scale_on_overbudget = false;
    } else if (arg == "--max-scale-promotions") {
      config.max_objective_scale_promotions =
          parseInt(requireValue(arg), arg);
    } else if (arg == "--random-initial-alpha-radius") {
      config.randomize_initial_alphas = true;
      config.initial_alpha_random_radius = parseLong(requireValue(arg), arg);
    } else if (arg == "--random-initial-alpha-seed") {
      config.initial_alpha_random_seed = static_cast<unsigned int>(
          std::strtoul(requireValue(arg).c_str(), nullptr, 10));
    } else if (arg == "--bk-storage") {
      config.bk_storage = requireValue(arg);
    } else if (arg == "--bk-mmap-dir") {
      config.bk_mmap_dir = requireValue(arg);
    } else if (arg == "--bk-mmap-advise") {
      config.bk_mmap_advise = requireValue(arg);
    } else {
      throw std::runtime_error("unknown argument: " + arg);
    }
  }
  if (config.partition_count <= 0) {
    throw std::runtime_error("--partitions must be positive");
  }
  if (config.max_iterations <= 0) {
    throw std::runtime_error("--max-iterations must be positive");
  }
  if (config.schedule_levels <= 0) {
    throw std::runtime_error("--schedule-levels must be positive");
  }
  if (config.schedule_start <= 0) {
    throw std::runtime_error("--schedule-start must be positive");
  }
  if (config.objective_scale <= 0) {
    throw std::runtime_error("--objective-scale must be positive");
  }
  if (config.directed && config.symmetric_streaming) {
    throw std::runtime_error(
        "--directed and --symmetric-streaming are mutually exclusive");
  }
  return config;
}

void configureBkStorage(const Config &config) {
  if (!config.bk_storage.empty()) {
    setenv("MCPD3_BK_STORAGE", config.bk_storage.c_str(), /*overwrite=*/1);
  }
  if (!config.bk_mmap_dir.empty()) {
    setenv("MCPD3_BK_MMAP_DIR", config.bk_mmap_dir.c_str(), /*overwrite=*/1);
  }
  if (!config.bk_mmap_advise.empty()) {
    setenv("MCPD3_BK_MMAP_ADVISE", config.bk_mmap_advise.c_str(),
           /*overwrite=*/1);
  }
}

int scaleIntCapacity(int value, long factor, bool saturate) {
  if (factor == 1 || value == 0) {
    return value;
  }
  const long int_max = std::numeric_limits<int>::max();
  const long int_min = std::numeric_limits<int>::min();
  const long value_long = static_cast<long>(value);
  const bool int_overflow =
      (value > 0 && factor > int_max / value_long) ||
      (value < 0 && factor > int_min / value_long);
  if (int_overflow && saturate) {
    return value > 0 ? std::numeric_limits<int>::max()
                     : std::numeric_limits<int>::min();
  }
  if (int_overflow) {
    throw std::overflow_error("objective scale exceeds int range");
  }
  const long scaled = static_cast<long>(value) * factor;
  return static_cast<int>(scaled);
}

void scaleGraph(mcpd3::MinCutGraph *graph, long factor, bool saturate) {
  if (factor == 1) {
    return;
  }
  for (auto &capacity : graph->arc_capacities) {
    capacity = scaleIntCapacity(capacity, factor, saturate);
  }
  for (auto &capacity : graph->terminal_capacities) {
    capacity = scaleIntCapacity(capacity, factor, saturate);
  }
}

mcpd3::DualDecompositionOptions makeOptions(const Config &config) {
  mcpd3::DualDecompositionOptions options;
  options.num_optimization_scales = config.schedule_levels;
  options.max_iteration_count = config.max_iterations;
  options.initial_step_size = config.schedule_start;
  options.patience = config.patience;
  options.legacy_patience = config.legacy_patience;
  options.use_momentum = config.use_momentum;
  options.enable_group_stopping = config.enable_group_stopping;
  options.track_primal_upper_bound = config.track_primal_upper_bound;
  options.emit_partition_packages = config.emit_partition_packages;
  options.verbose = config.verbose;
  options.objective_scale = config.objective_scale;
  options.thread_count = config.thread_count;
  options.regularization_scheme = config.regularization_scheme;
  options.regularization_budget_limit = config.regularization_budget_limit;
  options.promote_objective_scale_on_overbudget =
      config.promote_objective_scale_on_overbudget;
  options.max_objective_scale_promotions =
      config.max_objective_scale_promotions;
  options.randomize_initial_alphas = config.randomize_initial_alphas;
  options.initial_alpha_random_radius = config.initial_alpha_random_radius;
  options.initial_alpha_random_seed = config.initial_alpha_random_seed;
  return options;
}

void printMemory(const std::string &prefix, const MemorySnapshot &snapshot) {
  std::cout << prefix << "_vm_kb " << snapshot.vm_kb << "\n";
  std::cout << prefix << "_rss_kb " << snapshot.rss_kb << "\n";
}

void printConfig(const Config &config) {
  std::cout << "benchmark native_monolith\n";
  std::cout << "dimacs_path " << config.dimacs_path << "\n";
  std::cout << "partition_count " << config.partition_count << "\n";
  std::cout << "schedule_start " << config.schedule_start << "\n";
  std::cout << "schedule_levels " << config.schedule_levels << "\n";
  std::cout << "max_iterations " << config.max_iterations << "\n";
  std::cout << "objective_scale " << config.objective_scale << "\n";
  std::cout << "thread_count " << config.thread_count << "\n";
  std::cout << "directed " << config.directed << "\n";
  std::cout << "symmetric_streaming " << config.symmetric_streaming << "\n";
  std::cout << "track_primal_upper_bound "
            << config.track_primal_upper_bound << "\n";
  std::cout << "emit_partition_packages " << config.emit_partition_packages
            << "\n";
  std::cout << "regularization "
            << regularizationSchemeName(config.regularization_scheme) << "\n";
  std::cout << "regularization_budget_limit "
            << config.regularization_budget_limit << "\n";
  std::cout << "promote_objective_scale_on_overbudget "
            << config.promote_objective_scale_on_overbudget << "\n";
  std::cout << "max_objective_scale_promotions "
            << config.max_objective_scale_promotions << "\n";
  std::cout << "randomize_initial_alphas "
            << config.randomize_initial_alphas << "\n";
  std::cout << "initial_alpha_random_radius "
            << config.initial_alpha_random_radius << "\n";
  std::cout << "initial_alpha_random_seed "
            << config.initial_alpha_random_seed << "\n";
  std::cout << "bk_storage "
            << (config.bk_storage.empty() ? "env_or_default"
                                          : config.bk_storage)
            << "\n";
  std::cout << "bk_mmap_dir "
            << (config.bk_mmap_dir.empty() ? "env_or_default"
                                           : config.bk_mmap_dir)
            << "\n";
  std::cout << "bk_mmap_advise "
            << (config.bk_mmap_advise.empty() ? "env_or_default"
                                              : config.bk_mmap_advise)
            << "\n";
}

void printFinal(const Timing &timing, const Config &config,
                const mcpd3::DualDecomposition &dual_decomp,
                double peak_rss_kb) {
  const bool agreed = dual_decomp.getLastDisagreementCount() == 0;
  std::cout << "status " << (agreed ? 0 : 1) << "\n";
  std::cout << "status_name " << (agreed ? "agreement" : "not_converged")
            << "\n";
  std::cout << "final_objective_raw "
            << dual_decomp.getBestCertifiedLowerBoundRaw() << "\n";
  std::cout << "final_objective "
            << dual_decomp.getBestCertifiedLowerBound() << "\n";
  std::cout << "best_certified_lower_bound_raw "
            << dual_decomp.getBestCertifiedLowerBoundRaw() << "\n";
  std::cout << "best_certified_lower_bound "
            << dual_decomp.getBestCertifiedLowerBound() << "\n";
  std::cout << "best_regularized_objective_raw "
            << dual_decomp.getBestRegularizedObjectiveRaw() << "\n";
  std::cout << "best_regularized_objective "
            << dual_decomp.getBestRegularizedObjective() << "\n";
  std::cout << "best_upper_bound_raw " << dual_decomp.getBestUpperBoundRaw()
            << "\n";
  std::cout << "best_upper_bound " << dual_decomp.getBestUpperBound() << "\n";
  std::cout << "current_upper_bound_raw "
            << dual_decomp.getCurrentUpperBoundRaw() << "\n";
  std::cout << "objective_scale " << dual_decomp.getScale() << "\n";
  std::cout << "objective_scale_promotions "
            << dual_decomp.getObjectiveScalePromotionCount() << "\n";
  std::cout << "total_iterations "
            << dual_decomp.getTotalOptimizationIterations() << "\n";
  std::cout << "final_disagreement_count "
            << dual_decomp.getLastDisagreementCount() << "\n";
  std::cout << "final_disagreement_norm_sq "
            << dual_decomp.getLastDisagreementNormSq() << "\n";
  std::cout << "final_regularization_budget "
            << dual_decomp.getLastRegularizationBudget() << "\n";
  std::cout << "final_regularization_contribution "
            << dual_decomp.getLastRegularizationContribution() << "\n";
  std::cout << "final_regularization_anchor_sink_count "
            << dual_decomp.getLastRegularizationAnchorSinkCount() << "\n";
  std::cout << "final_regularization_active_sink_count "
            << dual_decomp.getLastRegularizationActiveSinkCount() << "\n";
  std::cout << "timing_total_wall_us " << timing.total_wall_us << "\n";
  std::cout << "timing_read_graph_wall_us " << timing.read_graph_wall_us
            << "\n";
  std::cout << "timing_scale_graph_wall_us " << timing.scale_graph_wall_us
            << "\n";
  std::cout << "timing_construct_wall_us " << timing.construct_wall_us
            << "\n";
  std::cout << "timing_solve_wall_us " << timing.solve_wall_us << "\n";
  std::cout << "timing_solver_inner_solve_wall_us "
            << dual_decomp.getTotalSolveLoopTime() << "\n";
  std::cout << "memory_peak_observed_rss_kb " << peak_rss_kb << "\n";
  std::cout << "saturate_capacity_overflow "
            << config.saturate_capacity_overflow << "\n";
}

} // namespace

int main(int argc, char **argv) {
  try {
    std::cout << std::unitbuf;
    const auto config = parseArgs(argc, argv);
    configureBkStorage(config);
    printConfig(config);

    Timing timing;
    const auto total_start = std::chrono::steady_clock::now();
    auto peak_rss_kb = memorySnapshot().rss_kb;
    printMemory("memory_initial", memorySnapshot());

    mcpd3::MinCutGraph graph;
    const auto read_start = std::chrono::steady_clock::now();
    if (config.directed) {
      graph = mcpd3::read_dimacs_directed_streaming(config.dimacs_path);
    } else if (config.symmetric_streaming) {
      graph = mcpd3::read_dimacs_symmetric_streaming(config.dimacs_path);
    } else {
      graph = mcpd3::read_dimacs(config.dimacs_path);
    }
    timing.read_graph_wall_us = elapsedUs(read_start);
    const auto after_read_memory = memorySnapshot();
    peak_rss_kb = std::max(peak_rss_kb, after_read_memory.rss_kb);
    printMemory("memory_after_read_graph", after_read_memory);

    std::cout << "graph_node_count " << graph.nnode << "\n";
    std::cout << "graph_arc_count " << graph.narc << "\n";

    const auto scale_start = std::chrono::steady_clock::now();
    scaleGraph(&graph, config.objective_scale,
               config.saturate_capacity_overflow);
    timing.scale_graph_wall_us = elapsedUs(scale_start);
    const auto after_scale_memory = memorySnapshot();
    peak_rss_kb = std::max(peak_rss_kb, after_scale_memory.rss_kb);
    printMemory("memory_after_scale_graph", after_scale_memory);

    std::unique_ptr<mcpd3::DualDecomposition> dual_decomp;
    const auto construct_start = std::chrono::steady_clock::now();
    dual_decomp = std::make_unique<mcpd3::DualDecomposition>(
        config.partition_count, std::move(graph), makeOptions(config));
    timing.construct_wall_us = elapsedUs(construct_start);
    const auto after_construct_memory = memorySnapshot();
    peak_rss_kb = std::max(peak_rss_kb, after_construct_memory.rss_kb);
    printMemory("memory_after_construct", after_construct_memory);

    const auto solve_start = std::chrono::steady_clock::now();
    dual_decomp->solve();
    timing.solve_wall_us = elapsedUs(solve_start);
    const auto after_solve_memory = memorySnapshot();
    peak_rss_kb = std::max(peak_rss_kb, after_solve_memory.rss_kb);
    printMemory("memory_after_solve", after_solve_memory);

    timing.total_wall_us = elapsedUs(total_start);
    printFinal(timing, config, *dual_decomp, peak_rss_kb);
  } catch (const std::exception &e) {
    std::cerr << "mcpd3_native_monolith_benchmark failed: " << e.what()
              << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
