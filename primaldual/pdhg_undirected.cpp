#include <primaldual/pdhg_undirected.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace mcpd3 {
namespace {

using Clock = std::chrono::steady_clock;
using Decimal = boost::multiprecision::cpp_dec_float_100;

long double elapsedSeconds(Clock::time_point start) {
  return std::chrono::duration<long double>(Clock::now() - start).count();
}

long double integerToLongDouble(const Objective &value) {
#if defined(MCPD_CAPACITY_MODE_32)
  return static_cast<long double>(value);
#else
  const std::string text = integer_to_string(value);
  std::size_t consumed = 0;
  const long double converted = std::stold(text, &consumed);
  if (consumed != text.size() || !std::isfinite(converted)) {
    throw std::overflow_error("PDHG objective cannot be represented as long double");
  }
  return converted;
#endif
}

long double nonnegativeIntegerToLowerLongDouble(const Objective &value) {
  long double converted = integerToLongDouble(value);
#if !defined(MCPD_CAPACITY_MODE_32)
  const Decimal exact(integer_to_string(value));
  if (Decimal(converted) > exact) {
    converted = std::nextafter(converted, 0.0L);
  }
#endif
  return converted;
}

long double nonnegativeIntegerToUpperLongDouble(const Objective &value) {
  long double converted = integerToLongDouble(value);
#if !defined(MCPD_CAPACITY_MODE_32)
  const Decimal exact(integer_to_string(value));
  if (Decimal(converted) < exact) {
    converted = std::nextafter(converted,
                               std::numeric_limits<long double>::infinity());
  }
#endif
  return converted;
}

long double clamp(long double value, long double minimum,
                  long double maximum) {
  return std::max(minimum, std::min(value, maximum));
}

struct SweepResult {
  Objective value = 0;
  std::vector<std::uint8_t> source_side;
};

struct DualEvaluation {
  long double raw_lower_bound = 0.0L;
  long double safe_lower_bound = 0.0L;
  long double conservation_residual = 0.0L;
};

void validateOptions(const PdhgOptions &options) {
  if (!std::isfinite(options.tau) || !std::isfinite(options.sigma) ||
      !std::isfinite(options.theta) ||
      !std::isfinite(options.step_size_scale) ||
      !std::isfinite(options.step_balance) ||
      !std::isfinite(options.capacity_quantum) ||
      !std::isfinite(options.absolute_gap_tolerance) ||
      !std::isfinite(options.relative_gap_tolerance) ||
      !std::isfinite(options.time_limit_seconds) ||
      !std::isfinite(options.stagnation_tolerance) ||
      !std::isfinite(options.lower_bound_safety_factor)) {
    throw std::invalid_argument("PDHG options must be finite");
  }
  if (options.check_interval == 0) {
    throw std::invalid_argument("PDHG check interval must be positive");
  }
  if (options.tau < 0.0L || options.sigma < 0.0L ||
      options.step_size_scale <= 0.0L || options.step_balance <= 0.0L) {
    throw std::invalid_argument("PDHG step sizes must be nonnegative and the scale positive");
  }
  if (options.theta < 0.0L || options.theta > 1.0L) {
    throw std::invalid_argument("PDHG theta must be in [0,1]");
  }
  if (options.capacity_quantum < 0.0L ||
      options.absolute_gap_tolerance < 0.0L ||
      options.relative_gap_tolerance < 0.0L ||
      options.time_limit_seconds < 0.0L ||
      options.stagnation_tolerance < 0.0L ||
      options.lower_bound_safety_factor < 0.0L) {
    throw std::invalid_argument("PDHG tolerances and limits must be nonnegative");
  }
}

} // namespace

const char *pdhg_termination_reason_name(PdhgTerminationReason reason) {
  switch (reason) {
  case PdhgTerminationReason::ExactCertificate:
    return "exact_certificate";
  case PdhgTerminationReason::ApproximateGap:
    return "approximate_gap";
  case PdhgTerminationReason::IterationLimit:
    return "iteration_limit";
  case PdhgTerminationReason::TimeLimit:
    return "time_limit";
  case PdhgTerminationReason::Stagnation:
    return "stagnation";
  }
  return "unknown";
}

PdhgUndirectedMinCutSolver::PdhgUndirectedMinCutSolver(
    const MinCutGraph &graph) {
  if (graph.nnode < 0 || graph.narc < 0 ||
      graph.arcs.size() != 2U * static_cast<std::size_t>(graph.narc) ||
      graph.arc_capacities.size() !=
          2U * static_cast<std::size_t>(graph.narc) ||
      graph.terminal_capacities.size() !=
          static_cast<std::size_t>(graph.nnode)) {
    throw std::invalid_argument("malformed MinCutGraph storage");
  }

  node_count_ = graph.nnode;
  source_ = node_count_;
  sink_ = node_count_ + 1;
  const int total_vertices = node_count_ + 2;
  std::vector<int> degree(static_cast<std::size_t>(total_vertices), 0);

  auto add_edge = [&](int source, int target, const Objective &capacity) {
    if (capacity < 0) {
      throw std::invalid_argument("PDHG requires nonnegative capacities");
    }
    if (capacity == 0 || source == target) {
      return;
    }
    edge_sources_.push_back(source);
    edge_targets_.push_back(target);
    edge_exact_capacities_.push_back(capacity);
    edge_capacities_.push_back(nonnegativeIntegerToLowerLongDouble(capacity));
    ++degree[static_cast<std::size_t>(source)];
    ++degree[static_cast<std::size_t>(target)];
  };

  edge_sources_.reserve(static_cast<std::size_t>(graph.narc + graph.nnode));
  edge_targets_.reserve(edge_sources_.capacity());
  edge_exact_capacities_.reserve(edge_sources_.capacity());
  edge_capacities_.reserve(edge_sources_.capacity());

  for (int edge = 0; edge < graph.narc; ++edge) {
    const int source = graph.arcs[2 * edge];
    const int target = graph.arcs[2 * edge + 1];
    if (source < 0 || source >= node_count_ || target < 0 ||
        target >= node_count_) {
      throw std::invalid_argument("PDHG edge endpoint is outside the graph");
    }
    const Capacity forward = graph.arc_capacities[2 * edge];
    const Capacity reverse = graph.arc_capacities[2 * edge + 1];
    if (forward < 0 || reverse < 0) {
      throw std::invalid_argument("PDHG requires nonnegative arc capacities");
    }
    if (forward != reverse) {
      throw std::invalid_argument(
          "PDHG undirected solver requires equal forward/reverse capacities");
    }
    add_edge(source, target, widen_capacity(forward));
  }

  for (int node = 0; node < node_count_; ++node) {
    const Capacity terminal = graph.terminal_capacities[node];
    if (terminal > 0) {
      add_edge(source_, node, widen_capacity(terminal));
    } else if (terminal < 0) {
      add_edge(node, sink_, checked_subtract(Objective{0},
                                             widen_capacity(terminal),
                                             "PDHG terminal magnitude overflow"));
    }
  }

  maximum_degree_ =
      node_count_ == 0
          ? 0
          : *std::max_element(degree.begin(), degree.begin() + node_count_);
  adjacency_offsets_.assign(static_cast<std::size_t>(total_vertices + 1), 0);
  for (const int vertex : edge_sources_) {
    ++adjacency_offsets_[static_cast<std::size_t>(vertex + 1)];
  }
  for (const int vertex : edge_targets_) {
    ++adjacency_offsets_[static_cast<std::size_t>(vertex + 1)];
  }
  std::partial_sum(adjacency_offsets_.begin(), adjacency_offsets_.end(),
                   adjacency_offsets_.begin());
  adjacency_edges_.resize(2 * edge_sources_.size());
  std::vector<std::size_t> cursor = adjacency_offsets_;
  for (std::size_t edge = 0; edge < edge_sources_.size(); ++edge) {
    adjacency_edges_[cursor[static_cast<std::size_t>(edge_sources_[edge])]++] =
        edge;
    adjacency_edges_[cursor[static_cast<std::size_t>(edge_targets_[edge])]++] =
        edge;
  }
}

PdhgResult PdhgUndirectedMinCutSolver::solve(
    const PdhgOptions &options) const {
  validateOptions(options);
  const auto start = Clock::now();
  const std::size_t total_vertices = static_cast<std::size_t>(node_count_ + 2);
  const std::size_t edge_count = edge_sources_.size();
  const long double automatic_step =
      maximum_degree_ == 0
          ? options.step_size_scale
          : options.step_size_scale /
                std::sqrt(2.0L * static_cast<long double>(maximum_degree_));
  const long double tau =
      options.tau > 0.0L ? options.tau : automatic_step / options.step_balance;
  const long double sigma =
      options.sigma > 0.0L ? options.sigma
                           : automatic_step * options.step_balance;

  std::vector<long double> x(total_vertices, 0.5L);
  std::vector<long double> x_new(total_vertices, 0.5L);
  std::vector<long double> x_bar(total_vertices, 0.5L);
  std::vector<long double> divergence(total_vertices, 0.0L);
  std::vector<long double> lower_divergence(total_vertices, 0.0L);
  std::vector<long double> p(edge_count, 0.0L);
  std::vector<long double> x_average;
  std::vector<long double> p_average;
  if (options.use_ergodic_primal) {
    x_average = x;
  }
  if (options.use_ergodic_dual) {
    p_average.assign(edge_count, 0.0L);
  }
  x[source_] = x_new[source_] = x_bar[source_] = 1.0L;
  x[sink_] = x_new[sink_] = x_bar[sink_] = 0.0L;
  if (!x_average.empty()) {
    x_average[source_] = 1.0L;
    x_average[sink_] = 0.0L;
  }

  std::vector<std::size_t> sweep_order(static_cast<std::size_t>(node_count_));
  std::iota(sweep_order.begin(), sweep_order.end(), 0);

  auto sweep = [&](const std::vector<long double> &labels) {
    std::sort(sweep_order.begin(), sweep_order.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                if (labels[lhs] != labels[rhs]) {
                  return labels[lhs] > labels[rhs];
                }
                return lhs < rhs;
              });
    std::vector<std::uint8_t> side(total_vertices, 0);
    side[static_cast<std::size_t>(source_)] = 1;
    Objective current = 0;
    for (std::size_t edge = adjacency_offsets_[source_];
         edge < adjacency_offsets_[source_ + 1]; ++edge) {
      current = checked_add(
          current, edge_exact_capacities_[adjacency_edges_[edge]],
          "PDHG sweep objective overflow");
    }
    SweepResult best{current, std::vector<std::uint8_t>(
                                  static_cast<std::size_t>(node_count_), 0)};
    for (const std::size_t node_index : sweep_order) {
      const int node = static_cast<int>(node_index);
      for (std::size_t offset = adjacency_offsets_[node];
           offset < adjacency_offsets_[node + 1]; ++offset) {
        const std::size_t edge = adjacency_edges_[offset];
        const int other = edge_sources_[edge] == node ? edge_targets_[edge]
                                                       : edge_sources_[edge];
        if (side[static_cast<std::size_t>(other)] != 0) {
          current = checked_subtract(current, edge_exact_capacities_[edge],
                                     "PDHG sweep objective underflow");
        } else {
          current = checked_add(current, edge_exact_capacities_[edge],
                                "PDHG sweep objective overflow");
        }
      }
      side[node_index] = 1;
      if (current < best.value) {
        best.value = current;
        std::copy(side.begin(), side.begin() + node_count_,
                  best.source_side.begin());
      }
    }
    return best;
  };

  auto compute_divergence = [&](const std::vector<long double> &dual,
                                std::vector<long double> *output) {
    std::fill(output->begin(), output->end(), 0.0L);
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
      (*output)[static_cast<std::size_t>(edge_sources_[edge])] += dual[edge];
      (*output)[static_cast<std::size_t>(edge_targets_[edge])] -= dual[edge];
    }
  };

  auto evaluate_dual = [&](const std::vector<long double> &dual) {
    compute_divergence(dual, &divergence);
    DualEvaluation evaluation;
    evaluation.raw_lower_bound = divergence[source_];
    long double residual_squared = 0.0L;
    for (int node = 0; node < node_count_; ++node) {
      evaluation.raw_lower_bound += std::min(0.0L, divergence[node]);
      residual_squared += divergence[node] * divergence[node];
    }
    evaluation.conservation_residual = std::sqrt(residual_squared);

    std::fill(lower_divergence.begin(), lower_divergence.end(), 0.0L);
    long double dual_magnitude = 0.0L;
    const long double negative_infinity =
        -std::numeric_limits<long double>::infinity();
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
      const std::size_t source =
          static_cast<std::size_t>(edge_sources_[edge]);
      const std::size_t target =
          static_cast<std::size_t>(edge_targets_[edge]);
      lower_divergence[source] = std::nextafter(
          lower_divergence[source] + dual[edge], negative_infinity);
      lower_divergence[target] = std::nextafter(
          lower_divergence[target] - dual[edge], negative_infinity);
      dual_magnitude += std::fabs(dual[edge]);
    }
    long double safe = lower_divergence[source_];
    for (int node = 0; node < node_count_; ++node) {
      safe = std::nextafter(safe + std::min(0.0L, lower_divergence[node]),
                            negative_infinity);
    }
    const long double margin =
        options.lower_bound_safety_factor *
        std::numeric_limits<long double>::epsilon() *
        std::max(1.0L, dual_magnitude);
    evaluation.safe_lower_bound =
        std::nextafter(safe - margin, negative_infinity);
    return evaluation;
  };

  auto fractional_objective = [&](const std::vector<long double> &labels) {
    long double value = 0.0L;
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
      value += edge_capacities_[edge] *
               std::fabs(labels[static_cast<std::size_t>(edge_sources_[edge])] -
                         labels[static_cast<std::size_t>(edge_targets_[edge])]);
    }
    return value;
  };

  PdhgResult result;
  result.effective_tau = tau;
  result.effective_sigma = sigma;
  bool has_upper_bound = false;
  bool has_lower_bound = false;
  std::size_t stagnant_checks = 0;
  Objective previous_best_cut = 0;
  long double last_x_change = 0.0L;
  long double last_p_change = 0.0L;

  auto check = [&](std::size_t iteration) {
    const SweepResult current_sweep = sweep(x);
    SweepResult check_best = current_sweep;
    if (!x_average.empty()) {
      SweepResult average_sweep = sweep(x_average);
      if (average_sweep.value < check_best.value) {
        check_best = std::move(average_sweep);
      }
    }
    const bool cut_improved =
        !has_upper_bound || check_best.value < result.best_cut_value;
    if (cut_improved) {
      has_upper_bound = true;
      result.best_cut_value = check_best.value;
      result.source_side = std::move(check_best.source_side);
      result.best_cut_first_iteration = iteration;
    }

    DualEvaluation dual = evaluate_dual(p);
    if (!p_average.empty()) {
      const DualEvaluation average_dual = evaluate_dual(p_average);
      if (average_dual.safe_lower_bound > dual.safe_lower_bound) {
        dual = average_dual;
      }
    }
    const bool lower_improved =
        !has_lower_bound ||
        dual.safe_lower_bound >
            result.safe_lower_bound + options.stagnation_tolerance;
    if (!has_lower_bound ||
        dual.safe_lower_bound > result.safe_lower_bound) {
      has_lower_bound = true;
      result.safe_lower_bound = dual.safe_lower_bound;
      result.best_dual_lower_bound = dual.raw_lower_bound;
    }

    const long double upper =
        nonnegativeIntegerToUpperLongDouble(result.best_cut_value);
    result.certified_gap = upper - result.safe_lower_bound;
    if (result.certified_gap < 0.0L) {
      throw std::runtime_error("PDHG safe lower bound exceeds its cut upper bound");
    }
    const long double relative_gap =
        result.certified_gap / std::max(1.0L, std::fabs(upper));
    PdhgIterationStats stats;
    stats.iteration = iteration;
    stats.elapsed_seconds = elapsedSeconds(start);
    stats.fractional_objective = fractional_objective(x);
    stats.current_cut_value = current_sweep.value;
    stats.best_cut_value = result.best_cut_value;
    stats.current_dual_lower_bound = dual.raw_lower_bound;
    stats.best_dual_lower_bound = result.best_dual_lower_bound;
    stats.safe_lower_bound = result.safe_lower_bound;
    stats.certified_gap = result.certified_gap;
    stats.relative_gap = relative_gap;
    stats.flow_conservation_residual = dual.conservation_residual;
    stats.x_change = last_x_change;
    stats.p_change = last_p_change;
    if (options.record_history) {
      result.history.push_back(stats);
    }
    if (options.verbose) {
      std::cout << "pdhg iteration=" << iteration
                << " elapsed_seconds=" << static_cast<double>(stats.elapsed_seconds)
                << " fractional_objective="
                << static_cast<double>(stats.fractional_objective)
                << " current_cut=" << integer_to_string(stats.current_cut_value)
                << " best_cut=" << integer_to_string(stats.best_cut_value)
                << " current_lower="
                << static_cast<double>(stats.current_dual_lower_bound)
                << " safe_lower=" << static_cast<double>(stats.safe_lower_bound)
                << " gap=" << static_cast<double>(stats.certified_gap)
                << " relative_gap=" << static_cast<double>(stats.relative_gap)
                << " conservation_residual="
                << static_cast<double>(stats.flow_conservation_residual)
                << " x_change=" << static_cast<double>(stats.x_change)
                << " p_change=" << static_cast<double>(stats.p_change) << '\n';
    }

    if (options.capacity_quantum > 0.0L &&
        result.certified_gap < options.capacity_quantum) {
      result.certified_exact = true;
      result.termination_reason = PdhgTerminationReason::ExactCertificate;
      return true;
    }
    const long double approximate_tolerance =
        options.absolute_gap_tolerance +
        options.relative_gap_tolerance * std::max(1.0L, std::fabs(upper));
    if (result.certified_gap <= approximate_tolerance) {
      result.termination_reason = PdhgTerminationReason::ApproximateGap;
      return true;
    }
    if (options.time_limit_seconds > 0.0L &&
        stats.elapsed_seconds >= options.time_limit_seconds) {
      result.termination_reason = PdhgTerminationReason::TimeLimit;
      return true;
    }
    if (options.stagnation_checks > 0) {
      const bool upper_improved =
          iteration == 0 || result.best_cut_value < previous_best_cut;
      if (upper_improved || lower_improved) {
        stagnant_checks = 0;
      } else {
        ++stagnant_checks;
      }
      previous_best_cut = result.best_cut_value;
      if (stagnant_checks >= options.stagnation_checks) {
        result.termination_reason = PdhgTerminationReason::Stagnation;
        return true;
      }
    }
    return false;
  };

  bool stopped = check(0);
  std::size_t iteration = 0;
  for (; !stopped && iteration < options.max_iterations; ++iteration) {
    long double p_change_squared = 0.0L;
    for (std::size_t edge = 0; edge < edge_count; ++edge) {
      const long double previous = p[edge];
      const long double updated =
          previous +
          sigma *
              (x_bar[static_cast<std::size_t>(edge_sources_[edge])] -
               x_bar[static_cast<std::size_t>(edge_targets_[edge])]);
      p[edge] = clamp(updated, -edge_capacities_[edge], edge_capacities_[edge]);
      const long double change = p[edge] - previous;
      p_change_squared += change * change;
    }

    compute_divergence(p, &divergence);
    long double x_change_squared = 0.0L;
    for (int node = 0; node < node_count_; ++node) {
      x_new[node] = clamp(x[node] - tau * divergence[node], 0.0L, 1.0L);
      const long double change = x_new[node] - x[node];
      x_change_squared += change * change;
      x_bar[node] = x_new[node] + options.theta * change;
    }
    x_new[source_] = x_bar[source_] = 1.0L;
    x_new[sink_] = x_bar[sink_] = 0.0L;
    x.swap(x_new);
    last_x_change = std::sqrt(x_change_squared);
    last_p_change = std::sqrt(p_change_squared);

    const long double average_count = static_cast<long double>(iteration + 1);
    if (!x_average.empty()) {
      for (std::size_t vertex = 0; vertex < total_vertices; ++vertex) {
        x_average[vertex] += (x[vertex] - x_average[vertex]) / average_count;
      }
    }
    if (!p_average.empty()) {
      for (std::size_t edge = 0; edge < edge_count; ++edge) {
        p_average[edge] += (p[edge] - p_average[edge]) / average_count;
      }
    }

    const std::size_t completed = iteration + 1;
    if (completed % options.check_interval == 0 ||
        completed == options.max_iterations) {
      stopped = check(completed);
    }
  }
  result.iterations = iteration;
  result.elapsed_seconds = elapsedSeconds(start);
  if (!stopped) {
    result.termination_reason = PdhgTerminationReason::IterationLimit;
  }
  return result;
}

} // namespace mcpd3
