#pragma once

#include <graph/mcgraph.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace mcpd3 {

enum class PdhgTerminationReason {
  ExactCertificate,
  ApproximateGap,
  IterationLimit,
  TimeLimit,
  Stagnation,
};

const char *pdhg_termination_reason_name(PdhgTerminationReason reason);

struct PdhgOptions {
  std::size_t max_iterations = 100000;
  std::size_t check_interval = 100;
  long double tau = 0.0L;
  long double sigma = 0.0L;
  long double theta = 1.0L;
  long double step_size_scale = 0.99L;
  long double step_balance = 1.0L;
  long double capacity_quantum = 1.0L;
  long double absolute_gap_tolerance = 1e-6L;
  long double relative_gap_tolerance = 1e-6L;
  long double time_limit_seconds = 0.0L;
  std::size_t stagnation_checks = 0;
  long double stagnation_tolerance = 0.0L;
  long double lower_bound_safety_factor = 64.0L;
  bool use_ergodic_primal = false;
  bool use_ergodic_dual = false;
  bool record_history = false;
  bool verbose = false;
};

struct PdhgIterationStats {
  std::size_t iteration = 0;
  long double elapsed_seconds = 0.0L;
  long double fractional_objective = 0.0L;
  Objective current_cut_value = 0;
  Objective best_cut_value = 0;
  long double current_dual_lower_bound = 0.0L;
  long double best_dual_lower_bound = 0.0L;
  long double safe_lower_bound = 0.0L;
  long double certified_gap = 0.0L;
  long double relative_gap = 0.0L;
  long double flow_conservation_residual = 0.0L;
  long double x_change = 0.0L;
  long double p_change = 0.0L;
};

struct PdhgResult {
  std::vector<std::uint8_t> source_side;
  Objective best_cut_value = 0;
  long double best_dual_lower_bound = 0.0L;
  long double safe_lower_bound = 0.0L;
  long double certified_gap = 0.0L;
  std::size_t iterations = 0;
  std::size_t best_cut_first_iteration = 0;
  long double elapsed_seconds = 0.0L;
  long double effective_tau = 0.0L;
  long double effective_sigma = 0.0L;
  bool certified_exact = false;
  PdhgTerminationReason termination_reason =
      PdhgTerminationReason::IterationLimit;
  std::vector<PdhgIterationStats> history;
};

class PdhgUndirectedMinCutSolver {
public:
  explicit PdhgUndirectedMinCutSolver(const MinCutGraph &graph);

  PdhgResult solve(const PdhgOptions &options = {}) const;

  int node_count() const { return node_count_; }
  std::size_t edge_count() const { return edge_sources_.size(); }
  int maximum_degree() const { return maximum_degree_; }

private:
  int node_count_ = 0;
  int source_ = 0;
  int sink_ = 0;
  int maximum_degree_ = 0;
  std::vector<int> edge_sources_;
  std::vector<int> edge_targets_;
  std::vector<long double> edge_capacities_;
  std::vector<Objective> edge_exact_capacities_;
  std::vector<std::size_t> adjacency_offsets_;
  std::vector<std::size_t> adjacency_edges_;
};

} // namespace mcpd3
