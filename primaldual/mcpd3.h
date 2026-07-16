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

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <decomp/constraint.h>
#include <graph/mcgraph.h>
#include <maxflow/graph.h>
#include <measure/timer.h>

namespace mcpd3 {

enum class CanonicalCutSelection {
  SOLVER_DEFAULT,
  MINIMUM_LABELS,
  MAXIMUM_LABELS
};

enum class ReferenceCutSelection {
  CLOSEST_EXACT,
  EXACT_REFERENCE_IF_OPTIMAL
};

inline bool primaldual_timing_enabled() {
  const char *value = std::getenv("MCPD3_SOLVER_TIMING");
  return value != nullptr && value[0] != '\0' && value[0] != '0';
}

class PrimalDualMinCutSolver {
public:
  using MaxflowGraph =
      Graph</*captype=*/Capacity, /*tcaptype=*/TerminalResidual,
            /*flowtype=*/Objective>;

  struct WarmState {
    std::vector<Capacity> v_flow;
    std::vector<NodeFlow> d_flow;
    std::vector<int> x;
    bool is_first_iteration = true;
    bool is_first_iteration_of_new_scale = true;
    bool has_solution = false;
    Objective mincut_value = 0;
    std::vector<Lagrange> cached_lagrange_multipliers;
    std::vector<Lagrange> cached_last_lagrange_multipliers;
    Capacity regularization_str = 0;
    Objective last_regularization_budget = 0;
    Objective last_regularization_contribution = 0;
    long last_regularization_anchor_sink_count = 0;
    long last_regularization_active_sink_count = 0;
    std::vector<Objective> regularization_weights;
    MaxflowGraph::ReusableState maxflow_graph_state;
  };

  struct FlowWarmStart {
    std::vector<int> arcs;
    std::vector<Capacity> arc_capacities;
    std::vector<Capacity> terminal_capacities;
    std::vector<Capacity> v_flow;
    std::vector<NodeFlow> d_flow;
    std::vector<int> x;
  };

  PrimalDualMinCutSolver(int nnode, int narc, std::vector<int> &&arcs,
                         std::vector<Capacity> arc_capacities,
                         std::vector<Capacity> terminal_capacities)
      : nnode_(nnode), narc_(narc), arcs_(std::move(arcs)),
        arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)), v_flow_(narc_, 0),
        d_flow_(nnode_, 0), x_(nnode_, 0),
        incremental_changed_node_flags_(nnode_, 0),
        maxflow_graph_(nnode_, narc_),
        is_first_iteration_(true), is_first_iteration_of_new_scale_(true),
        has_solution_(false),
        canonical_cut_selection_(CanonicalCutSelection::SOLVER_DEFAULT),
        reference_cut_selection_(ReferenceCutSelection::CLOSEST_EXACT),
        force_full_mincut_recompute_(false),
        maxflow_changed_list_(128),
        regularization_str_(0),
        last_regularization_budget_(0), last_regularization_contribution_(0),
        last_regularization_anchor_sink_count_(0),
        last_regularization_active_sink_count_(0) {
    initializeMaxflowGraph();
  }

  template <typename InputCapacity,
            std::enable_if_t<!std::is_same_v<InputCapacity, Capacity>, int> = 0>
  PrimalDualMinCutSolver(int nnode, int narc, std::vector<int> &&arcs,
                         const std::vector<InputCapacity> &arc_capacities,
                         const std::vector<InputCapacity> &terminal_capacities)
      : PrimalDualMinCutSolver(
            nnode, narc, std::move(arcs),
            capacity_vector_from(arc_capacities),
            capacity_vector_from(terminal_capacities)) {}

  PrimalDualMinCutSolver(MinCutGraph min_cut_graph)
      : PrimalDualMinCutSolver(min_cut_graph.nnode, min_cut_graph.narc,
                               std::move(min_cut_graph.arcs),
                               std::move(min_cut_graph.arc_capacities),
                               std::move(min_cut_graph.terminal_capacities)) {}

  void setTrackArcFlowUpdates(bool enabled) {
    if (enabled == track_arc_flow_updates_) {
      return;
    }
    track_arc_flow_updates_ = enabled;
    if (enabled) {
      arc_flow_update_counts_.assign(static_cast<size_t>(narc_), 0);
    } else {
      arc_flow_update_counts_.clear();
      arc_flow_update_counts_.shrink_to_fit();
    }
  }

  const std::vector<std::uint64_t> &getArcFlowUpdateCounts() const {
    return arc_flow_update_counts_;
  }

  void resetArcFlowUpdateCounts() {
    std::fill(arc_flow_update_counts_.begin(), arc_flow_update_counts_.end(),
              std::uint64_t{0});
  }

  void decodeNarrowBand(const std::list<int> seeds, int rad) {
    // TODO: this function needs to be cleaned up and rewritten to be much more
    // memory/runtime efficient

    std::vector<bool> visited(nnode_, false);
    std::vector<int> dist(nnode_, 0);
    std::list<int> q;

    int index = 0;
    for (const auto &index : seeds) {
      q.emplace_back(index);
      visited[index] = true;
    }
    MaxflowGraph::arc_id a;
    auto nodes = maxflow_graph_.get_nodes();

    while (!q.empty()) {
      int u = q.front();
      q.pop_front();
      if (dist[u] >= rad) {
        continue;
      }
      const auto &_u = nodes[u];
      for (a = _u.first; a; a = a->next) {
        auto v = a->head;
        auto iv = std::distance(nodes, v);
        if (!visited[iv]) {
          dist[iv] = dist[u] + 1;
          visited[iv] = true;
          q.emplace_back(iv);
        }
      }
    }

    size_t num_arcs_in_decoding = 0;
    size_t num_nodes_in_decoding = 0;

    // setup capacities
    a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      if (!visited[s] || !visited[t]) {
        maxflow_graph_.set_rcap(a, 0);
        a = maxflow_graph_.get_next_arc(a);
        maxflow_graph_.set_rcap(a, 0);
        a = maxflow_graph_.get_next_arc(a);
        continue;
      }
      num_arcs_in_decoding++;
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      maxflow_graph_.set_rcap(a, forward_capacity);
      a = maxflow_graph_.get_next_arc(a);
      maxflow_graph_.set_rcap(a, backward_capacity);
      a = maxflow_graph_.get_next_arc(a);
    }
    for (int i = 0; i < nnode_; ++i) {
      if (!visited[i]) {
        maxflow_graph_.set_trcap(i, 0);
        continue;
      }
      num_nodes_in_decoding++;
      if (dist[i] < rad) {
        auto terminal_capacity = terminal_capacities_[i];
        if (terminal_capacity > 0) {
          maxflow_graph_.add_tweights(i, terminal_capacity, 0);
        } else {
          maxflow_graph_.add_tweights(i, 0, -terminal_capacity);
        }
      } else {
        if (x_[i] == 0) {
          maxflow_graph_.set_trcap(i, std::numeric_limits<int>::max() / 2);
        } else {
          maxflow_graph_.set_trcap(i, -std::numeric_limits<int>::max() / 2);
        }
      }
    }

    auto maxflow_val = maxflow_graph_.maxflow();

    for (int i = 0; i < nnode_; ++i) {
      if (!visited[i]) {
        continue;
      }
      if (dist[i] < rad) {
        x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK;
      }
    }

    printf("/////////////////////////////////////////\n");
    printf("//////////  DECODING STATS //////////////\n");
    printf("// num_arcs in decoding :  %8lu    //\n", num_arcs_in_decoding);
    printf("// num_nodes in decoding : %8lu    //\n", num_nodes_in_decoding);
    printf("// maxflow val :           %s    //\n",
           integer_to_string(maxflow_val).c_str());
    printf("/////////////////////////////////////////\n");

    computeMinCutValueInitial(); // TODO: is this needed?
  }

  template <int scale> void scaleProblem() { scaleProblem(scale); }

  void scaleProblem(long scale, bool saturate_capacity_overflow = false) {
    if (scale <= 0) {
      throw std::runtime_error("problem scale factor must be positive");
    }
    for (int i = 0; i < nnode_; ++i) {
      terminal_capacities_[i] =
          checked_scale_capacity(terminal_capacities_[i], scale,
                                 saturate_capacity_overflow);
      auto &flow = d_flow_[i];
      flow = checked_scale(flow, scale, "node flow scale overflow");
    }
    for (int i = 0; i < narc_; ++i) {
      auto &forward_capacity = arc_capacities_[2 * i + 0];
      auto &backward_capacity = arc_capacities_[2 * i + 1];
      forward_capacity =
          checked_scale_capacity(forward_capacity, scale,
                                 saturate_capacity_overflow);
      backward_capacity =
          checked_scale_capacity(backward_capacity, scale,
                                 saturate_capacity_overflow);
      auto &flow = v_flow_[i];
      flow = checked_scale_capacity(flow, scale, saturate_capacity_overflow);
    }
    // The BK residual also contains flow induced by regularization. Scaling it
    // would therefore implement scale * (F + R), not scale * F + R. Preserve
    // the explicit arc flow above, then rebuild an exact residual network from
    // that warm flow on the next solve while leaving R unchanged.
    maxflow_graph_.reset();
    initializeMaxflowGraph();
    maxflow_changed_list_.Reset();
    incremental_mincut_nodes_.clear();
    incremental_arcs_.clear();
    mincut_value_ =
        checked_scale(mincut_value_, scale, "objective scale promotion overflow");
    // NOTE: after changing scale, the capacities from previous and this scale
    // are at completely different values, hence incremental update of
    // mincut_value_ will not work properly
    is_first_iteration_ = true;
    is_first_iteration_of_new_scale_ = true;
  }

  void setRegularizationStrength(const Capacity &str) {
    if (str < 0) {
      throw std::runtime_error("regularization strength must be non-negative");
    }
    regularization_str_ = str;
  }

  void setCanonicalCutSelection(CanonicalCutSelection selection) {
    canonical_cut_selection_ = selection;
  }

  CanonicalCutSelection getCanonicalCutSelection() const {
    return canonical_cut_selection_;
  }

  void setReferenceCutLabels(std::vector<int> labels) {
    if (labels.size() != static_cast<size_t>(nnode_)) {
      throw std::runtime_error(
          "reference cut label count must match the local node count");
    }
    for (const int label : labels) {
      if (label != 0 && label != 1) {
        throw std::runtime_error("reference cut labels must be binary");
      }
    }
    reference_cut_labels_ = std::move(labels);
  }

  void clearReferenceCutLabels() { reference_cut_labels_.clear(); }

  bool hasReferenceCutLabels() const {
    return !reference_cut_labels_.empty();
  }

  void setReferenceCutSelection(ReferenceCutSelection selection) {
    reference_cut_selection_ = selection;
  }

  void setReferenceCutCheckInterval(long interval) {
    if (interval <= 0) {
      throw std::runtime_error(
          "reference cut check interval must be positive");
    }
    reference_cut_check_interval_ = interval;
  }

  long getReferenceDecodeCount() const { return reference_decode_count_; }
  long getReferenceCurrentCutHitCount() const {
    return reference_current_cut_hit_count_;
  }
  long getReferenceExactHitCount() const {
    return reference_exact_hit_count_;
  }
  long getReferenceClosureCount() const { return reference_closure_count_; }
  long getReferenceDecodeTimeMicroseconds() const {
    return reference_decode_time_us_;
  }

  void setForceFullMinCutRecompute(bool enabled) {
    force_full_mincut_recompute_ = enabled;
  }

  Objective maxflow() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      maxflow_graph_.set_rcap(a, forward_capacity);
      a = maxflow_graph_.get_next_arc(a);
      maxflow_graph_.set_rcap(a, backward_capacity);
      a = maxflow_graph_.get_next_arc(a);
    }
    for (int i = 0; i < nnode_; ++i) {
      auto terminal_capacity = terminal_capacities_[i];
      if (terminal_capacity > 0) {
        maxflow_graph_.add_tweights(i, terminal_capacity, 0);
      } else {
        maxflow_graph_.add_tweights(i, 0, -terminal_capacity);
      }
    }
    auto maxflow = maxflow_graph_.maxflow();
    maxflow_graph_.reset();
    initializeMaxflowGraph();
    return maxflow;
  }

  void solve() {
    if (is_first_iteration_) {
      auto init_time = time_lambda([&] {
        shrinkToFitDualDecompositionConstraints(); // memory optimization
        initializeFlow(); // finds a flow satisfying arc based lagrange
                          // multiplier complementary slackness conditions
      });
      if (primaldual_timing_enabled()) {
        printf("init_time: %ldms\n", init_time.count());
      }
    }
    cacheLagrangeMultipliers(); // optimization
    resetRegularizationDiagnostics();
    updateRegularizationAnchorsFromCurrentSolution();
    updateNodePotentials();     // finds which node based lagrange multiplier
                                // complementary slackness conditions are
                                // violated and sets source and sink capacities
                                // accordingly
    computeMaxflow();           // compute maxflow
    updateFlow();               // get updated flow
    updateMinCut();             // get updated min cut solution
    has_solution_ = true;

    // set flag indicating that incremental methods should be used hereafter
    if (is_first_iteration_ || is_first_iteration_of_new_scale_) {
      computeMinCutValueInitial(); // initialize min cut value to compute
                                   // incremental changes later
      is_first_iteration_of_new_scale_ = false;
      is_first_iteration_ = false;
    }
  }

  Objective getMaxFlowValue() const {
    Objective maxflow = 0;
    std::vector<NodeFlow> node_balance(nnode_, 0);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      const Capacity flow = v_flow_[i];
      node_balance[s] = checked_add(
          node_balance[s], node_flow_from_capacity(flow),
          "maxflow balance overflow");
      node_balance[t] = checked_subtract(
          node_balance[t], node_flow_from_capacity(flow),
          "maxflow balance overflow");
    }
    for (int i = 0; i < nnode_; ++i) {
      if (node_balance[i] > 0) {
        maxflow = checked_add(maxflow,
                              static_cast<Objective>(node_balance[i]),
                              "maxflow objective overflow");
      }
      // TODO: the imbalance needs to be accounted for
    }
    return maxflow;
  }

  Objective getMinCutValue() const { return mincut_value_; }
  Capacity getRegularizationStrength() const { return regularization_str_; }
  Objective getLastRegularizationBudget() const {
    return last_regularization_budget_;
  }
  Objective getLastRegularizationContribution() const {
    return last_regularization_contribution_;
  }
  long getLastRegularizationAnchorSinkCount() const {
    return last_regularization_anchor_sink_count_;
  }
  long getLastRegularizationActiveSinkCount() const {
    return last_regularization_active_sink_count_;
  }

  void addSourceDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_source;
    auto find_iter = std::find(dual_decomposition_local_indices_.begin(),
                               dual_decomposition_local_indices_.end(), index);
    if (find_iter != dual_decomposition_local_indices_.end()) {
      auto &constraint = dual_decomposition_constraints_[std::distance(
          dual_decomposition_local_indices_.begin(), find_iter)];
      constraint.source_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.source_arc_references.emplace_back(arc_reference);
      dual_decomposition_local_indices_.emplace_back(index);
      dual_decomposition_constraints_.emplace_back(constraint);
    }
  }

  void addTargetDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_target;
    auto find_iter = std::find(dual_decomposition_local_indices_.begin(),
                               dual_decomposition_local_indices_.end(), index);
    if (find_iter != dual_decomposition_local_indices_.end()) {
      auto &constraint = dual_decomposition_constraints_[std::distance(
          dual_decomposition_local_indices_.begin(), find_iter)];
      constraint.target_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.target_arc_references.emplace_back(arc_reference);
      dual_decomposition_local_indices_.emplace_back(index);
      dual_decomposition_constraints_.emplace_back(constraint);
    }
  }

  int getMinCutSolution(int index) const { return x_[index]; }

  struct MemoryEstimate {
    std::size_t bk_node_bytes = 0;
    std::size_t bk_arc_bytes = 0;
    std::size_t bk_total_bytes = 0;
    std::size_t solver_vector_bytes = 0;
    std::size_t total_bytes = 0;
  };

  static MemoryEstimate estimateMemoryBytes(int nnode, int narc) {
    using EstimateGraph =
        Graph</*captype=*/Capacity, /*tcaptype=*/TerminalResidual,
              /*flowtype=*/Objective>;
    MemoryEstimate estimate;
    estimate.bk_node_bytes = EstimateGraph::estimated_node_array_bytes(nnode);
    estimate.bk_arc_bytes = EstimateGraph::estimated_arc_array_bytes(narc);
    estimate.bk_total_bytes =
        estimate.bk_node_bytes + estimate.bk_arc_bytes;
    const auto arc_index_count = 2 * static_cast<std::size_t>(narc);
    const auto incremental_arc_index_count =
        static_cast<std::size_t>(narc);
    const auto arc_capacity_count = 3 * static_cast<std::size_t>(narc);
    const auto arc_change_flag_count = static_cast<std::size_t>(narc);
    const auto node_capacity_count = static_cast<std::size_t>(nnode);
    const auto node_flow_count = static_cast<std::size_t>(nnode);
    const auto node_label_count = static_cast<std::size_t>(nnode);
    const auto node_change_flag_count = static_cast<std::size_t>(nnode);
    const auto incremental_node_index_count =
        static_cast<std::size_t>(nnode);
    estimate.solver_vector_bytes =
        (arc_index_count + incremental_arc_index_count +
         incremental_node_index_count) * sizeof(int) +
        (arc_capacity_count + node_capacity_count) * sizeof(Capacity) +
        node_flow_count * sizeof(NodeFlow) +
        node_label_count * sizeof(int) +
        (node_change_flag_count + arc_change_flag_count) *
            sizeof(unsigned char);
    estimate.total_bytes =
        estimate.bk_total_bytes + estimate.solver_vector_bytes;
    return estimate;
  }

  void setMinCutSolution(const std::vector<bool> &new_solution) {
    std::copy(new_solution.begin(), new_solution.end(), x_.begin());
    computeMinCutValueInitial();
    has_solution_ = true;
  }

  void setMinCutSolution(const std::vector<int> &new_solution) {
    std::copy(new_solution.begin(), new_solution.end(), x_.begin());
    computeMinCutValueInitial();
    has_solution_ = true;
  }

  FlowWarmStart captureFlowWarmStart() const {
    if (!has_solution_) {
      throw std::runtime_error("cannot capture flow warm start before solve");
    }
    if (regularization_str_ != 0 ||
        std::any_of(regularization_weights_.begin(),
                    regularization_weights_.end(),
                    [](const Objective &weight) { return weight != 0; })) {
      throw std::runtime_error(
          "cannot capture flow warm start from a regularized solve");
    }
    return FlowWarmStart{arcs_, arc_capacities_, terminal_capacities_,
                         v_flow_, d_flow_, x_};
  }

  void restoreFlowWarmStart(const FlowWarmStart &state) {
    if (!is_first_iteration_ || has_solution_) {
      throw std::runtime_error(
          "flow warm start must be restored before the first solve");
    }
    if (regularization_str_ != 0) {
      throw std::runtime_error(
          "flow warm start requires regularization to be disabled");
    }
    if (state.arcs != arcs_) {
      throw std::runtime_error("flow warm start graph topology mismatch");
    }
    if (state.arc_capacities.size() != arc_capacities_.size() ||
        state.terminal_capacities.size() != terminal_capacities_.size() ||
        state.v_flow.size() != v_flow_.size() ||
        state.d_flow.size() != d_flow_.size() || state.x.size() != x_.size()) {
      throw std::runtime_error("flow warm start shape mismatch");
    }
    for (size_t i = 0; i < arc_capacities_.size(); ++i) {
      if (arc_capacities_[i] < state.arc_capacities[i]) {
        throw std::runtime_error(
            "flow warm start arc capacity decreased");
      }
    }
    for (size_t i = 0; i < terminal_capacities_.size(); ++i) {
      const Capacity old_capacity = state.terminal_capacities[i];
      const Capacity new_capacity = terminal_capacities_[i];
      if (old_capacity == 0) {
        continue;
      }
      const bool same_sign = (old_capacity > 0) == (new_capacity > 0);
      const Objective old_magnitude = absolute_capacity(old_capacity);
      const Objective new_magnitude = absolute_capacity(new_capacity);
      if (!same_sign || new_magnitude < old_magnitude) {
        throw std::runtime_error(
            "flow warm start terminal capacity is not monotone");
      }
    }

    v_flow_ = state.v_flow;
    d_flow_ = state.d_flow;
    x_ = state.x;
    has_solution_ = true;
  }

  void validateProblemCapacityReplacement(
      const std::vector<Capacity> &arc_capacities,
      const std::vector<Capacity> &terminal_capacities,
      bool preserve_flow_state = true,
      const Objective &flow_scale_numerator = 1,
      const Objective &flow_scale_denominator = 1) const {
    if (arc_capacities.size() != arc_capacities_.size()) {
      throw std::runtime_error("replacement arc capacity count mismatch");
    }
    if (terminal_capacities.size() != terminal_capacities_.size()) {
      throw std::runtime_error(
          "replacement terminal capacity count mismatch");
    }
    for (const Capacity &capacity : arc_capacities) {
      if (capacity < 0) {
        throw std::runtime_error(
            "replacement arc capacities must be non-negative");
      }
    }

    if (flow_scale_numerator <= 0 || flow_scale_denominator <= 0) {
      throw std::invalid_argument("flow scale ratio must be positive");
    }

    if (preserve_flow_state &&
        flow_scale_numerator != flow_scale_denominator) {
      for (size_t i = 0; i < arc_capacities.size(); ++i) {
        if (!capacities_have_ratio(
                arc_capacities_[i], arc_capacities[i],
                flow_scale_numerator, flow_scale_denominator)) {
          throw std::runtime_error(
              "flow scaling requires proportional arc capacities");
        }
      }
    }
  }

  void replaceProblemCapacities(const std::vector<Capacity> &arc_capacities,
                                const std::vector<Capacity> &terminal_capacities,
                                bool preserve_flow_state = true,
                                const Objective &flow_scale_numerator = 1,
                                const Objective &flow_scale_denominator = 1) {
    validateProblemCapacityReplacement(
        arc_capacities, terminal_capacities, preserve_flow_state,
        flow_scale_numerator, flow_scale_denominator);

    std::vector<Capacity> replacement_flow;
    if (preserve_flow_state &&
        flow_scale_numerator != flow_scale_denominator) {
      replacement_flow.reserve(v_flow_.size());
      for (const Capacity &flow : v_flow_) {
        replacement_flow.push_back(checked_scale_capacity_ratio(
            flow, flow_scale_numerator, flow_scale_denominator));
      }
    }

    arc_capacities_ = arc_capacities;
    terminal_capacities_ = terminal_capacities;
    if (!preserve_flow_state) {
      std::fill(v_flow_.begin(), v_flow_.end(), 0);
    } else if (!replacement_flow.empty()) {
      v_flow_ = std::move(replacement_flow);
    }
    std::fill(d_flow_.begin(), d_flow_.end(), 0);
    for (int i = 0; i < narc_; ++i) {
      const Capacity lower = -arc_capacities_[2 * i + 1];
      const Capacity upper = arc_capacities_[2 * i];
      v_flow_[i] = std::max(lower, std::min(upper, v_flow_[i]));
      const int source = arcs_[2 * i];
      const int target = arcs_[2 * i + 1];
      d_flow_[source] = checked_add(
          d_flow_[source], node_flow_from_capacity(v_flow_[i]),
          "warm-start node balance overflow");
      d_flow_[target] = checked_subtract(
          d_flow_[target], node_flow_from_capacity(v_flow_[i]),
          "warm-start node balance overflow");
    }

    maxflow_graph_.reset();
    initializeMaxflowGraph();
    maxflow_changed_list_.Reset();
    incremental_mincut_nodes_.clear();
    incremental_arcs_.clear();
    is_first_iteration_ = true;
    is_first_iteration_of_new_scale_ = true;
    mincut_value_ = 0;
    resetRegularizationDiagnostics();
  }

  template <typename InputCapacity,
            std::enable_if_t<!std::is_same_v<InputCapacity, Capacity>, int> = 0>
  void replaceProblemCapacities(
      const std::vector<InputCapacity> &arc_capacities,
      const std::vector<InputCapacity> &terminal_capacities,
      bool preserve_flow_state = true,
      const Objective &flow_scale_numerator = 1,
      const Objective &flow_scale_denominator = 1) {
    replaceProblemCapacities(capacity_vector_from(arc_capacities),
                             capacity_vector_from(terminal_capacities),
                             preserve_flow_state, flow_scale_numerator,
                             flow_scale_denominator);
  }

  WarmState captureWarmState() const {
    WarmState state;
    state.v_flow = v_flow_;
    state.d_flow = d_flow_;
    state.x = x_;
    state.is_first_iteration = is_first_iteration_;
    state.is_first_iteration_of_new_scale = is_first_iteration_of_new_scale_;
    state.has_solution = has_solution_;
    state.mincut_value = mincut_value_;
    state.cached_lagrange_multipliers = cached_lagrange_multipliers_;
    state.cached_last_lagrange_multipliers =
        cached_last_lagrange_multipliers_;
    state.regularization_str = regularization_str_;
    state.last_regularization_budget = last_regularization_budget_;
    state.last_regularization_contribution =
        last_regularization_contribution_;
    state.last_regularization_anchor_sink_count =
        last_regularization_anchor_sink_count_;
    state.last_regularization_active_sink_count =
        last_regularization_active_sink_count_;
    state.regularization_weights = regularization_weights_;
    state.maxflow_graph_state = maxflow_graph_.captureReusableState();
    return state;
  }

  void restoreWarmState(const WarmState &state) {
    if (state.v_flow.size() != static_cast<size_t>(narc_) ||
        state.d_flow.size() != static_cast<size_t>(nnode_) ||
        state.x.size() != static_cast<size_t>(nnode_)) {
      throw std::runtime_error("solver warm state shape does not match graph");
    }
    v_flow_ = state.v_flow;
    d_flow_ = state.d_flow;
    x_ = state.x;
    is_first_iteration_ = state.is_first_iteration;
    is_first_iteration_of_new_scale_ = state.is_first_iteration_of_new_scale;
    has_solution_ = state.has_solution;
    mincut_value_ = state.mincut_value;
    cached_lagrange_multipliers_ = state.cached_lagrange_multipliers;
    cached_last_lagrange_multipliers_ =
        state.cached_last_lagrange_multipliers;
    regularization_str_ = state.regularization_str;
    last_regularization_budget_ = state.last_regularization_budget;
    last_regularization_contribution_ =
        state.last_regularization_contribution;
    last_regularization_anchor_sink_count_ =
        state.last_regularization_anchor_sink_count;
    last_regularization_active_sink_count_ =
        state.last_regularization_active_sink_count;
    regularization_weights_ = state.regularization_weights;
    incremental_mincut_nodes_.clear();
    incremental_arcs_.clear();
    maxflow_changed_list_.Reset();
    dual_decomposition_local_indices_set_.clear();
    for (const auto &index : dual_decomposition_local_indices_) {
      dual_decomposition_local_indices_set_.insert(index);
    }
    maxflow_graph_.restoreReusableState(state.maxflow_graph_state);
  }

private:

  void resetRegularizationDiagnostics() {
    last_regularization_budget_ = 0;
    last_regularization_contribution_ = 0;
    last_regularization_anchor_sink_count_ = 0;
    last_regularization_active_sink_count_ = 0;
    if (regularization_weights_.size() !=
        dual_decomposition_local_indices_.size()) {
      regularization_weights_.assign(dual_decomposition_local_indices_.size(),
                                     Objective{0});
    }
  }

  Lagrange lagrangeMultiplierTerm(size_t constraint_index) const {
    const auto &constraint = dual_decomposition_constraints_[constraint_index];
    Lagrange lagrange_multiplier_term = 0;
    for (const auto &arc_reference : constraint.source_arc_references) {
      lagrange_multiplier_term = checked_subtract(
          lagrange_multiplier_term, arc_reference->alpha,
          "lagrange multiplier term overflow");
    }
    for (const auto &arc_reference : constraint.target_arc_references) {
      lagrange_multiplier_term = checked_add(
          lagrange_multiplier_term, arc_reference->alpha,
          "lagrange multiplier term overflow");
    }
    return lagrange_multiplier_term;
  }

  TerminalResidual regularizationTerm(size_t constraint_index) const {
    if (regularization_weights_.empty()) {
      return 0;
    }
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
    return terminal_residual_from_capacity(
        capacity_from_integer(regularization_weights_[constraint_index]));
#else
    return static_cast<TerminalResidual>(
        regularization_weights_[constraint_index]);
#endif
  }

  Objective updateRegularizationAnchorsFromCurrentSolution() {
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
    if (regularization_str_ <= 0) {
      std::fill(regularization_weights_.begin(), regularization_weights_.end(),
                Objective{0});
    } else {
      for (Objective &weight : regularization_weights_) {
        if (weight != 0) {
          weight = widen_capacity(regularization_str_);
        }
      }
    }
#endif
    if (has_solution_) {
      for (size_t i = 0; i < dual_decomposition_local_indices_.size(); ++i) {
        if (cached_lagrange_multipliers_[i] ==
            cached_last_lagrange_multipliers_[i]) {
          continue;
        }
        const int local_index = dual_decomposition_local_indices_[i];
#if defined(MCPD_LEGACY_32BIT_DD_REPLAY)
        regularization_weights_[i] =
            x_[local_index] ? widen_capacity(regularization_str_)
                            : Objective{0};
#else
        if (regularization_str_ > 0 && x_[local_index]) {
          regularization_weights_[i] = checked_add(
              regularization_weights_[i], widen_capacity(regularization_str_),
              "cumulative scaled epsilon regularization overflow");
        }
#endif
      }
    }
    Objective budget = 0;
    for (size_t i = 0; i < dual_decomposition_local_indices_.size(); ++i) {
      if (regularization_weights_[i] > 0) {
        budget = checked_add(budget, regularization_weights_[i],
                             "scaled epsilon regularization budget overflow");
        last_regularization_anchor_sink_count_++;
      }
    }
    last_regularization_budget_ = budget;
    return budget;
  }

  void computeMinCutValueInitial() {
    mincut_value_ = 0;
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      if (x_[s] == 0 && x_[t] == 1) {
        mincut_value_ = checked_add(
            mincut_value_, widen_capacity(forward_capacity),
            "mincut objective overflow");
      } else if (x_[s] == 1 && x_[t] == 0) {
        mincut_value_ = checked_add(
            mincut_value_, widen_capacity(backward_capacity),
            "mincut objective overflow");
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      auto terminal_capacity = terminal_capacities_[i];
      if (x_[i] == 0 && terminal_capacity < 0) {
        mincut_value_ = checked_add(
            mincut_value_, absolute_capacity(terminal_capacity),
            "mincut terminal objective overflow");
      } else if (x_[i] == 1 && terminal_capacity > 0) {
        mincut_value_ = checked_add(
            mincut_value_, widen_capacity(terminal_capacity),
            "mincut terminal objective overflow");
      }
    }
    // add dual decomposition node potential terms (when/if applicable)
    size_t i = 0;
    for (const auto &constraint : dual_decomposition_constraints_) {
      (void)constraint;
      const Lagrange lagrange_multiplier_term = lagrangeMultiplierTerm(i);
      if (x_[dual_decomposition_local_indices_[i]]) {
        mincut_value_ = checked_add(
            mincut_value_, widen_lagrange(lagrange_multiplier_term),
            "mincut lagrange objective overflow");
      }
      ++i;
    }
  }

  void initializeMaxflowGraph() {
    maxflow_graph_.add_node(nnode_);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      maxflow_graph_.add_edge(s, t, 0, 0);
    }
  }

  std::pair<Capacity, Capacity>
  arcGradients(const Capacity &forward_capacity,
               const Capacity &backward_capacity,
               const Capacity &flow) const {
    Capacity pos = checked_add(flow, forward_capacity,
                               "forward residual capacity overflow");
    Capacity neg = checked_subtract(flow, backward_capacity,
                                    "backward residual capacity overflow");
    return {pos, neg};
  }

  TerminalResidual nodeGradient(const Capacity &terminal_capacity,
                                const NodeFlow &flow) const {
    return checked_add(flow, terminal_residual_from_capacity(terminal_capacity),
                       "terminal residual capacity overflow");
  }

  void initializeFlow() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      auto flow = v_flow_[i];
      auto [pos, neg] = arcGradients(forward_capacity, backward_capacity, flow);
      Capacity new_flow = 0;
      Capacity dfp, dfn;
      if (pos < 0 || neg > 0) {
        dfp = std::min(pos, Capacity(0));
        dfn = std::max(neg, Capacity(0));
        if (absolute_capacity(dfn) > absolute_capacity(dfp)) {
          new_flow = -dfn;
        } else {
          new_flow = -dfp;
        }
      }
      if (new_flow != 0) {
        v_flow_[i] = checked_add(v_flow_[i], new_flow, "arc flow overflow");
        d_flow_[s] = checked_add(d_flow_[s], node_flow_from_capacity(new_flow),
                                 "node flow balance overflow");
        d_flow_[t] = checked_subtract(d_flow_[t], node_flow_from_capacity(new_flow),
                                      "node flow balance overflow");
      }
      std::tie(pos, neg) =
          arcGradients(forward_capacity, backward_capacity, v_flow_[i]);
      assert(pos >= 0 && neg <= 0);
      maxflow_graph_.set_rcap(a, pos);
      a = maxflow_graph_.get_next_arc(a);
      maxflow_graph_.set_rcap(a, -neg);
      a = maxflow_graph_.get_next_arc(a);
    }
  }

  void updateNodeTerminal(int i, TerminalResidual pos, bool do_update) {
    if (do_update) {
      auto existing_pos = maxflow_graph_.get_trcap(i);
      pos = checked_add(pos, existing_pos,
                        "terminal residual capacity overflow");
      if (!is_first_iteration_) {
        if (existing_pos != pos) {
          maxflow_graph_.set_trcap(i, pos);
          maxflow_graph_.mark_node(i);
        }
      } else if (is_first_iteration_) {
        maxflow_graph_.set_trcap(i, pos);
      }
    } else {
      if (!is_first_iteration_) {
        auto stored_pos = maxflow_graph_.get_trcap(i);
        if (stored_pos != pos) {
          maxflow_graph_.set_trcap(i, pos);
          maxflow_graph_.mark_node(i);
        }
      } else if (is_first_iteration_) {
        maxflow_graph_.set_trcap(i, pos);
      }
    }
  }

  void updateNodePotentials() {
    if (is_first_iteration_) {
      updateNodePotentialsInitial();
    } else {
      updateNodePotentialsIncremental();
    }
  }

  void updateNodePotentialsInitial() {
    // add mincut node potential terms
    for (int i = 0; i < nnode_; ++i) {
      auto pos = nodeGradient(terminal_capacities_[i], d_flow_[i]);
      updateNodeTerminal(i, pos, false);
    }
    updateDualDecompositionNodePotentials();
  }

  void updateNodePotentialsIncremental() {
    // add mincut node potential terms
    for (const int i : incremental_mincut_nodes_) {
      if ( dual_decomposition_local_indices_set_.find(i) != dual_decomposition_local_indices_set_.end() ){
        continue; // will be handled separately below
      }
      auto pos = nodeGradient(terminal_capacities_[i], d_flow_[i]);
      updateNodeTerminal(i, pos, false);
    }
    for (size_t cache_index = 0;
         cache_index < dual_decomposition_local_indices_.size();
         ++cache_index) {
      const int i = dual_decomposition_local_indices_[cache_index];
      auto pos = nodeGradient(terminal_capacities_[i], d_flow_[i]);
      pos = checked_add(pos, cached_lagrange_multipliers_[cache_index],
                        "terminal lagrange capacity overflow");
      pos = checked_add(pos, regularizationTerm(cache_index),
                        "terminal regularization capacity overflow");
      updateNodeTerminal(i, pos, false);
    }
  }

  void updateDualDecompositionNodePotentials() {
    // add dual decomposition node potential terms (when/if applicable)
    size_t cache_index = 0;
    for (const auto &index : dual_decomposition_local_indices_) {
      auto pos = cached_lagrange_multipliers_[cache_index];
      pos = checked_add(pos, regularizationTerm(cache_index),
                        "terminal regularization capacity overflow");
      cache_index++;
      updateNodeTerminal(index, pos, true);
    }
  }

  void updateMinCut() {
    const bool check_reference =
        !reference_cut_labels_.empty() &&
        reference_cut_schedule_count_++ % reference_cut_check_interval_ == 0;
    if (check_reference) {
      const auto decode_time = time_lambda([&] {
        ++reference_decode_count_;
        updateMinCutInitial();
        if (x_ == reference_cut_labels_) {
          ++reference_current_cut_hit_count_;
        } else if (isReferenceCutOptimal()) {
          ++reference_exact_hit_count_;
          x_ = reference_cut_labels_;
        } else if (reference_cut_selection_ ==
                   ReferenceCutSelection::CLOSEST_EXACT) {
          ++reference_closure_count_;
          updateReferenceGuidedMinCut();
        }
        computeMinCutValueInitial();
      });
      reference_decode_time_us_ += decode_time.count();
    } else if (canonical_cut_selection_ !=
               CanonicalCutSelection::SOLVER_DEFAULT) {
      updateCanonicalMinCut();
      computeMinCutValueInitial();
    } else if (is_first_iteration_ || force_full_mincut_recompute_) {
      updateMinCutInitial();
      if (!is_first_iteration_) {
        computeMinCutValueInitial();
      }
    } else {
      updateMinCutIncremental();
    }
    updateRegularizationContribution();
  }

  void updateCanonicalMinCut() {
    std::vector<unsigned char> reached(static_cast<size_t>(nnode_), 0);
    std::deque<int> queue;
    auto enqueue = [&](int node) {
      if (!reached[static_cast<size_t>(node)]) {
        reached[static_cast<size_t>(node)] = 1;
        queue.push_back(node);
      }
    };

    if (canonical_cut_selection_ == CanonicalCutSelection::MAXIMUM_LABELS) {
      for (int node = 0; node < nnode_; ++node) {
        if (maxflow_graph_.get_trcap(node) > 0) {
          enqueue(node);
        }
      }
    } else {
      for (int node = 0; node < nnode_; ++node) {
        if (maxflow_graph_.get_trcap(node) < 0) {
          enqueue(node);
        }
      }
    }

    auto nodes = maxflow_graph_.get_nodes();
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop_front();
      for (MaxflowGraph::arc_id arc = nodes[node].first; arc;
           arc = arc->next) {
        const bool traversable =
            canonical_cut_selection_ == CanonicalCutSelection::MAXIMUM_LABELS
                ? maxflow_graph_.get_rcap(arc) > 0
                : maxflow_graph_.get_rcap(arc->sister) > 0;
        if (traversable) {
          enqueue(static_cast<int>(std::distance(nodes, arc->head)));
        }
      }
    }

    for (int node = 0; node < nnode_; ++node) {
      if (canonical_cut_selection_ == CanonicalCutSelection::MAXIMUM_LABELS) {
        // Source-reachable nodes form the minimum source-side min-cut.
        x_[node] = reached[static_cast<size_t>(node)] ? 0 : 1;
      } else {
        // Its dual: nodes that can reach the sink must remain sink-side.
        x_[node] = reached[static_cast<size_t>(node)] ? 1 : 0;
      }
    }
  }

  bool isReferenceCutOptimal() {
    for (int node = 0; node < nnode_; ++node) {
      const Objective terminal = maxflow_graph_.get_trcap(node);
      const int label = reference_cut_labels_[static_cast<size_t>(node)];
      if ((terminal > 0 && label != 0) || (terminal < 0 && label != 1)) {
        return false;
      }
    }

    auto arc = maxflow_graph_.get_first_arc();
    for (int index = 0; index < maxflow_graph_.get_arc_num(); ++index) {
      MaxflowGraph::node_id source = -1;
      MaxflowGraph::node_id target = -1;
      maxflow_graph_.get_arc_ends(arc, source, target);
      if (maxflow_graph_.get_rcap(arc) > 0 &&
          reference_cut_labels_[static_cast<size_t>(source)] == 0 &&
          reference_cut_labels_[static_cast<size_t>(target)] == 1) {
        return false;
      }
      arc = maxflow_graph_.get_next_arc(arc);
    }
    return true;
  }

  void updateReferenceGuidedMinCut() {
    using ArcId = MaxflowGraph::arc_id;
    struct DfsFrame {
      int node = -1;
      ArcId next = nullptr;
    };

    auto nodes = maxflow_graph_.get_nodes();
    std::vector<unsigned char> visited(static_cast<size_t>(nnode_), 0);
    std::vector<int> finish_order;
    finish_order.reserve(static_cast<size_t>(nnode_));
    std::vector<DfsFrame> dfs;
    for (int root = 0; root < nnode_; ++root) {
      if (visited[static_cast<size_t>(root)]) {
        continue;
      }
      visited[static_cast<size_t>(root)] = 1;
      dfs.push_back(DfsFrame{root, nodes[root].first});
      while (!dfs.empty()) {
        auto &frame = dfs.back();
        ArcId arc = frame.next;
        while (arc != nullptr) {
          const int target =
              static_cast<int>(std::distance(nodes, arc->head));
          if (maxflow_graph_.get_rcap(arc) > 0 &&
              !visited[static_cast<size_t>(target)]) {
            break;
          }
          arc = arc->next;
        }
        if (arc == nullptr) {
          finish_order.push_back(frame.node);
          dfs.pop_back();
          continue;
        }
        frame.next = arc->next;
        const int target =
            static_cast<int>(std::distance(nodes, arc->head));
        visited[static_cast<size_t>(target)] = 1;
        dfs.push_back(DfsFrame{target, nodes[target].first});
      }
    }

    std::vector<int> component(static_cast<size_t>(nnode_), -1);
    std::vector<int> pending;
    int component_count = 0;
    for (auto iter = finish_order.rbegin(); iter != finish_order.rend();
         ++iter) {
      const int root = *iter;
      if (component[static_cast<size_t>(root)] != -1) {
        continue;
      }
      component[static_cast<size_t>(root)] = component_count;
      pending.push_back(root);
      while (!pending.empty()) {
        const int node = pending.back();
        pending.pop_back();
        for (ArcId arc = nodes[node].first; arc != nullptr; arc = arc->next) {
          if (maxflow_graph_.get_rcap(arc->sister) <= 0) {
            continue;
          }
          const int predecessor =
              static_cast<int>(std::distance(nodes, arc->head));
          if (component[static_cast<size_t>(predecessor)] == -1) {
            component[static_cast<size_t>(predecessor)] = component_count;
            pending.push_back(predecessor);
          }
        }
      }
      ++component_count;
    }

    if (nnode_ == std::numeric_limits<int>::max()) {
      throw std::overflow_error(
          "reference-guided cut is too large for closure capacities");
    }
    const Capacity implication_capacity = capacity_from_integer(nnode_ + 1);
    std::vector<Objective> component_terminal(
        static_cast<size_t>(component_count), 0);
    std::vector<unsigned char> forced_source(
        static_cast<size_t>(component_count), 0);
    std::vector<unsigned char> forced_sink(
        static_cast<size_t>(component_count), 0);
    for (int node = 0; node < nnode_; ++node) {
      const int id = component[static_cast<size_t>(node)];
      component_terminal[static_cast<size_t>(id)] = checked_add(
          component_terminal[static_cast<size_t>(id)],
          reference_cut_labels_[static_cast<size_t>(node)] == 0 ? Objective{1}
                                                                : Objective{-1},
          "reference closure terminal capacity overflow");
      const Objective terminal = maxflow_graph_.get_trcap(node);
      if (terminal > 0) {
        forced_source[static_cast<size_t>(id)] = 1;
      } else if (terminal < 0) {
        forced_sink[static_cast<size_t>(id)] = 1;
      }
    }
    for (int id = 0; id < component_count; ++id) {
      if (forced_source[static_cast<size_t>(id)] &&
          forced_sink[static_cast<size_t>(id)]) {
        throw std::runtime_error(
            "residual component is forced to both terminals");
      }
      if (forced_source[static_cast<size_t>(id)]) {
        component_terminal[static_cast<size_t>(id)] = checked_add(
            component_terminal[static_cast<size_t>(id)],
            widen_capacity(implication_capacity),
            "reference closure terminal capacity overflow");
      }
      if (forced_sink[static_cast<size_t>(id)]) {
        component_terminal[static_cast<size_t>(id)] = checked_subtract(
            component_terminal[static_cast<size_t>(id)],
            widen_capacity(implication_capacity),
            "reference closure terminal capacity overflow");
      }
    }

    std::vector<std::pair<int, int>> implications;
    implications.reserve(static_cast<size_t>(narc_) * 2);
    ArcId arc = maxflow_graph_.get_first_arc();
    for (int edge = 0; edge < narc_; ++edge) {
      const int source = arcs_[2 * edge];
      const int target = arcs_[2 * edge + 1];
      const int source_component = component[static_cast<size_t>(source)];
      const int target_component = component[static_cast<size_t>(target)];
      if (source_component != target_component &&
          maxflow_graph_.get_rcap(arc) > 0) {
        implications.emplace_back(source_component, target_component);
      }
      arc = maxflow_graph_.get_next_arc(arc);
      if (source_component != target_component &&
          maxflow_graph_.get_rcap(arc) > 0) {
        implications.emplace_back(target_component, source_component);
      }
      arc = maxflow_graph_.get_next_arc(arc);
    }

    if (implications.size() >
        static_cast<size_t>(std::numeric_limits<int>::max())) {
      throw std::overflow_error(
          "reference-guided closure has too many implications");
    }
    MaxflowGraph closure_graph(component_count,
                               static_cast<int>(implications.size()));
    closure_graph.add_node(component_count);
    for (const auto &[source, target] : implications) {
      closure_graph.add_edge(source, target, implication_capacity, 0);
    }
    for (int id = 0; id < component_count; ++id) {
      const Objective terminal = component_terminal[static_cast<size_t>(id)];
      if (terminal > 0) {
        closure_graph.add_tweights(id, terminal, 0);
      } else if (terminal < 0) {
        closure_graph.add_tweights(id, 0, -terminal);
      }
    }
    (void)closure_graph.maxflow();
    for (int node = 0; node < nnode_; ++node) {
      x_[node] =
          closure_graph.what_segment(component[static_cast<size_t>(node)]) ==
                  MaxflowGraph::SINK
              ? 1
              : 0;
    }
  }

  void updateRegularizationContribution() {
    last_regularization_contribution_ = 0;
    last_regularization_active_sink_count_ = 0;
    if (regularization_weights_.empty()) {
      return;
    }
    for (size_t i = 0; i < dual_decomposition_local_indices_.size(); ++i) {
      const int local_index = dual_decomposition_local_indices_[i];
      if (regularization_weights_[i] > 0 && x_[local_index]) {
        last_regularization_contribution_ = checked_add(
            last_regularization_contribution_, regularization_weights_[i],
            "regularization contribution overflow");
        last_regularization_active_sink_count_++;
      }
    }
  }

  void updateMinCutInitial() {
    for (int i = 0; i < nnode_; ++i) {
      x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
    }
  }

  void updateMinCutIncremental() {
    // update dual decomposition node potential terms (when/if applicable)
    size_t cache_index = 0;
    for (const auto &index : dual_decomposition_local_indices_) {
      auto x_i_new =
          maxflow_graph_.what_segment(index) == MaxflowGraph::SINK ? 1 : 0;
      Lagrange lagrange_multiplier_term =
          cached_lagrange_multipliers_[cache_index];
      Lagrange last_lagrange_multiplier_term =
          cached_last_lagrange_multipliers_[cache_index];
      if (last_lagrange_multiplier_term == lagrange_multiplier_term &&
          x_i_new == x_[index]) {
        cache_index++;
        continue;
      } else {
        cache_index++;
      }
      if (x_[index]) {
        mincut_value_ = checked_subtract(
            mincut_value_, widen_lagrange(last_lagrange_multiplier_term),
            "incremental mincut lagrange overflow");
      }
      if (x_i_new) {
        mincut_value_ = checked_add(
            mincut_value_, widen_lagrange(lagrange_multiplier_term),
            "incremental mincut lagrange overflow");
      }
    }

    // Mark every node whose cut label actually changed. Keeping the old labels
    // intact until all incident edges are evaluated lets one endpoint own an
    // edge changed at both ends, without allocating a per-solve hash set.
    for (const int i : incremental_mincut_nodes_) {
      if (incremental_changed_node_flags_[static_cast<size_t>(i)] != 0) {
        continue;
      }
      const int x_i_new =
          maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
      if (x_i_new != x_[i]) {
        incremental_changed_node_flags_[static_cast<size_t>(i)] = 1;
      }
    }

    // Update node and arc terms that may have changed.
    auto nodes = maxflow_graph_.get_nodes();
    MaxflowGraph::arc_id first_arc = maxflow_graph_.get_first_arc();
    for (const int i : incremental_mincut_nodes_) {
      auto &changed =
          incremental_changed_node_flags_[static_cast<size_t>(i)];
      if (changed != 1) {
        continue;
      }
      const int x_i_new =
          maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;

      // process terminals
      auto terminal_capacity = terminal_capacities_[i];
      if (x_[i] == 0 && x_i_new == 1) {
        mincut_value_ = checked_add(
            mincut_value_, widen_capacity(terminal_capacity),
            "incremental mincut terminal overflow");
      }
      if (x_[i] == 1 && x_i_new == 0) {
        mincut_value_ = checked_subtract(
            mincut_value_, widen_capacity(terminal_capacity),
            "incremental mincut terminal overflow");
      }

      // processes each possible arc
      MaxflowGraph::arc_id a;
      const auto &node_i = nodes[i];
      for (a = node_i.first; a; a = a->next) {
        auto arc_index = std::distance(first_arc, a) / 2;
        int s = arcs_[2 * arc_index + 0];
        int t = arcs_[2 * arc_index + 1];
        if (s == t) {
          continue;
        }
        const bool source_changed =
            incremental_changed_node_flags_[static_cast<size_t>(s)] != 0;
        const bool target_changed =
            incremental_changed_node_flags_[static_cast<size_t>(t)] != 0;
        if (source_changed && target_changed && i != s) {
          continue;
        }
        auto forward_capacity = arc_capacities_[2 * arc_index + 0];
        auto backward_capacity = arc_capacities_[2 * arc_index + 1];
        auto x_s_new =
            maxflow_graph_.what_segment(s) == MaxflowGraph::SINK ? 1 : 0;
        auto x_t_new =
            maxflow_graph_.what_segment(t) == MaxflowGraph::SINK ? 1 : 0;
        if ((x_[s] == 0 && x_[t] == 1) && !(x_s_new == 0 && x_t_new == 1)) {
          mincut_value_ = checked_subtract(
              mincut_value_, widen_capacity(forward_capacity),
              "incremental mincut arc overflow");
        }
        if ((x_[s] == 1 && x_[t] == 0) && !(x_s_new == 1 && x_t_new == 0)) {
          mincut_value_ = checked_subtract(
              mincut_value_, widen_capacity(backward_capacity),
              "incremental mincut arc overflow");
        }
        if (!(x_[s] == 0 && x_[t] == 1) && (x_s_new == 0 && x_t_new == 1)) {
          mincut_value_ = checked_add(
              mincut_value_, widen_capacity(forward_capacity),
              "incremental mincut arc overflow");
        }
        if (!(x_[s] == 1 && x_[t] == 0) && (x_s_new == 1 && x_t_new == 0)) {
          mincut_value_ = checked_add(
              mincut_value_, widen_capacity(backward_capacity),
              "incremental mincut arc overflow");
        }
      }
      changed = 2;
    }

    for (const int i : incremental_mincut_nodes_) {
      auto &changed =
          incremental_changed_node_flags_[static_cast<size_t>(i)];
      if (changed == 0) {
        continue;
      }
      x_[i] =
          maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
      changed = 0;
    }
  }

  void updateFlowInitial() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      auto flow = v_flow_[i];
      auto [pos, neg] = arcGradients(forward_capacity, backward_capacity, flow);
      Capacity new_flow = checked_subtract(
          maxflow_graph_.get_rcap(a), pos, "arc flow delta overflow");
      recordArcFlowUpdate(i, new_flow);
      v_flow_[i] = checked_add(v_flow_[i], new_flow, "arc flow overflow");
      d_flow_[s] = checked_add(d_flow_[s], node_flow_from_capacity(new_flow),
                               "node flow balance overflow");
      d_flow_[t] = checked_subtract(
          d_flow_[t], node_flow_from_capacity(new_flow),
          "node flow balance overflow");
      a = maxflow_graph_.get_next_arc(a);
      a = maxflow_graph_.get_next_arc(a);
    }
  }

  void updateFlowIncremental() {
    MaxflowGraph::arc_id first_arc = maxflow_graph_.get_first_arc();
    for (const int i : incremental_arcs_) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      auto flow = v_flow_[i];
      auto [pos, neg] = arcGradients(forward_capacity, backward_capacity, flow);
      Capacity new_flow = checked_subtract(
          maxflow_graph_.get_rcap(first_arc + 2 * i), pos,
          "arc flow delta overflow");
      recordArcFlowUpdate(i, new_flow);
      v_flow_[i] = checked_add(v_flow_[i], new_flow, "arc flow overflow");
      d_flow_[s] = checked_add(d_flow_[s], node_flow_from_capacity(new_flow),
                               "node flow balance overflow");
      d_flow_[t] = checked_subtract(d_flow_[t], node_flow_from_capacity(new_flow),
                                    "node flow balance overflow");
    }
  }

  void updateFlow() {
    if (is_first_iteration_) {
      updateFlowInitial();
    } else {
      updateFlowIncremental();
    }
  }

  void recordArcFlowUpdate(int arc_index, const Capacity &flow_delta) {
    if (!track_arc_flow_updates_ || flow_delta == 0) {
      return;
    }
    auto &count = arc_flow_update_counts_[static_cast<size_t>(arc_index)];
    if (count == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error("arc flow update count overflow");
    }
    ++count;
  }

  void computeMaxflow() {
    if (is_first_iteration_) {
      maxflow_graph_.maxflow();
    } else {
      incremental_mincut_nodes_.clear();
      maxflow_graph_.maxflow(true, incremental_arcs_, &maxflow_changed_list_);

      // update incremental nodes
      MaxflowGraph::node_id *ptr;
      for (ptr = maxflow_changed_list_.ScanFirst(); ptr;
           ptr = maxflow_changed_list_.ScanNext()) {
        MaxflowGraph::node_id i = *ptr;
        maxflow_graph_.remove_from_changed_list(i);
        incremental_mincut_nodes_.emplace_back(i);
      }
      maxflow_changed_list_.Reset();
    }
  }

  void cacheLagrangeMultipliers() {
    if (is_first_iteration_) {
      cached_lagrange_multipliers_.resize(
          dual_decomposition_local_indices_.size());
      cached_last_lagrange_multipliers_.resize(
          dual_decomposition_local_indices_.size());
    }
    size_t cache_index = 0;
    for (const auto &constraint : dual_decomposition_constraints_) {
      Lagrange lagrange_multiplier_term = 0;
      Lagrange last_lagrange_multiplier_term = 0;
      for (const auto &arc_reference : constraint.source_arc_references) {
        lagrange_multiplier_term = checked_subtract(
            lagrange_multiplier_term, arc_reference->alpha,
            "lagrange multiplier cache overflow");
        last_lagrange_multiplier_term = checked_subtract(
            last_lagrange_multiplier_term, arc_reference->last_alpha,
            "lagrange multiplier cache overflow");
      }
      for (const auto &arc_reference : constraint.target_arc_references) {
        lagrange_multiplier_term = checked_add(
            lagrange_multiplier_term, arc_reference->alpha,
            "lagrange multiplier cache overflow");
        last_lagrange_multiplier_term = checked_add(
            last_lagrange_multiplier_term, arc_reference->last_alpha,
            "lagrange multiplier cache overflow");
      }
      cached_lagrange_multipliers_[cache_index] = lagrange_multiplier_term;
      cached_last_lagrange_multipliers_[cache_index] =
          last_lagrange_multiplier_term;
      cache_index++;
    }
  }

  void shrinkToFitDualDecompositionConstraints() {
    dual_decomposition_local_indices_.shrink_to_fit();
    for (auto &constraint : dual_decomposition_constraints_) {
      constraint.source_arc_references.shrink_to_fit();
      constraint.target_arc_references.shrink_to_fit();
    }
    for ( const auto &i : dual_decomposition_local_indices_ ) {
      dual_decomposition_local_indices_set_.insert(i);
    }
  }

  /**
   * data passed into solver
   */
  int nnode_;
  int narc_;
  std::vector<int> arcs_;
  std::vector<Capacity> arc_capacities_;
  std::vector<Capacity> terminal_capacities_;

  /**
   * data structures needed for solving primal dual problem
   */
  std::vector<Capacity> v_flow_; // flow on the arcs
  std::vector<std::uint64_t> arc_flow_update_counts_;
  bool track_arc_flow_updates_ = false;
  std::vector<NodeFlow> d_flow_; // flow balance on the nodes
  std::vector<int> x_;      // mincut solution
  std::vector<unsigned char> incremental_changed_node_flags_;
  MaxflowGraph maxflow_graph_; // graph used to compute maxflow
  bool is_first_iteration_;
  bool is_first_iteration_of_new_scale_;
  bool has_solution_;
  CanonicalCutSelection canonical_cut_selection_;
  ReferenceCutSelection reference_cut_selection_;
  std::vector<int> reference_cut_labels_;
  long reference_cut_check_interval_ = 1;
  long reference_cut_schedule_count_ = 0;
  long reference_decode_count_ = 0;
  long reference_current_cut_hit_count_ = 0;
  long reference_exact_hit_count_ = 0;
  long reference_closure_count_ = 0;
  long reference_decode_time_us_ = 0;
  bool force_full_mincut_recompute_;

  Block<MaxflowGraph::node_id> maxflow_changed_list_;
  std::vector<int> incremental_mincut_nodes_;
  std::vector<int> incremental_arcs_;
  Objective mincut_value_;

  /**
   * specific to dual decomposition
   */
  struct DualDecompositionConstraint {
    std::vector<DualDecompositionConstraintArcReference> source_arc_references;
    std::vector<DualDecompositionConstraintArcReference> target_arc_references;
  };
  std::vector<int> dual_decomposition_local_indices_;
  std::unordered_set<int> dual_decomposition_local_indices_set_;
  std::vector<DualDecompositionConstraint> dual_decomposition_constraints_;

  std::vector<Lagrange> cached_lagrange_multipliers_;
  std::vector<Lagrange> cached_last_lagrange_multipliers_;

  Capacity regularization_str_;
  Objective last_regularization_budget_ = 0;
  Objective last_regularization_contribution_ = 0;
  long last_regularization_anchor_sink_count_ = 0;
  long last_regularization_active_sink_count_ = 0;
  std::vector<Objective> regularization_weights_;
};

} // namespace mcpd3
