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
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <decomp/constraint.h>
#include <graph/mcgraph.h>
#include <maxflow/graph.h>

namespace mcpd3 {

class PrimalDualMinCutSolver {
public:
  PrimalDualMinCutSolver(int nnode, int narc, std::vector<int> &&arcs,
                         std::vector<int> arc_capacities,
                         std::vector<int> terminal_capacities)
      : nnode_(nnode), narc_(narc), arcs_(std::move(arcs)),
        arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)), v_flow_(narc_, 0),
        d_flow_(nnode_, 0), x_(nnode_, 0), maxflow_graph_(nnode_, narc_),
        is_first_iteration_(true), is_first_iteration_of_new_scale_(true),
        maxflow_changed_list_(128) {
    initializeMaxflowGraph();
  }

  PrimalDualMinCutSolver(MinCutGraph min_cut_graph)
      : PrimalDualMinCutSolver(min_cut_graph.nnode, min_cut_graph.narc,
                               std::move(min_cut_graph.arcs),
                               std::move(min_cut_graph.arc_capacities),
                               std::move(min_cut_graph.terminal_capacities)) {}

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
    printf("// maxflow val :           %8ld    //\n", maxflow_val);
    printf("/////////////////////////////////////////\n");

    computeMinCutValueInitial(); // TODO: is this needed?
  }

  template <int scale> void scaleProblem() {
    for (int i = 0; i < nnode_; ++i) {
      terminal_capacities_[i] *= scale;
      auto &flow = d_flow_[i];
      flow *= scale;
    }
    for (int i = 0; i < narc_; ++i) {
      auto &forward_capacity = arc_capacities_[2 * i + 0];
      auto &backward_capacity = arc_capacities_[2 * i + 1];
      forward_capacity *= scale;
      backward_capacity *= scale;
      auto &flow = v_flow_[i];
      flow *= scale;
    }
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int flow;
      flow = maxflow_graph_.get_rcap(a);
      maxflow_graph_.set_rcap(a, flow * scale);
      a = maxflow_graph_.get_next_arc(a);
      flow = maxflow_graph_.get_rcap(a);
      maxflow_graph_.set_rcap(a, flow * scale);
      a = maxflow_graph_.get_next_arc(a);
    }
    for (int i = 0; i < nnode_; ++i) {
      int flow;
      flow = maxflow_graph_.get_trcap(i);
      maxflow_graph_.set_trcap(i, flow * scale);
    }
    // NOTE: after changing scale, the capacities from previous and this scale
    // are at completely different values, hence incremental update of
    // mincut_value_ will not work properly
    is_first_iteration_of_new_scale_ = true;
  }

  long maxflow() {
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
      shrinkToFitDualDecompositionConstraints(); // memory optimization
      initializeFlow(); // finds a flow satisfying arc based lagrange multiplier
                        // complementary slackness conditions
    }
    cacheLagrangeMultipliers(); // optimization
    updateNodePotentials();     // finds which node based lagrange multiplier
                                // complementary slackness conditions are
                                // violated and sets source and sink capacities
                                // accordingly
    computeMaxflow();           // compute maxflow
    updateFlow();               // get updated flow
    updateMinCut();             // get updated min cut solution

    // set flag indicating that incremental methods should be used hereafter
    if (is_first_iteration_ || is_first_iteration_of_new_scale_) {
      computeMinCutValueInitial(); // initialize min cut value to compute
                                   // incremental changes later
      is_first_iteration_of_new_scale_ = false;
      is_first_iteration_ = false;
    }
  }

  long getMaxFlowValue() const {
    long maxflow = 0;
    std::vector<int> node_balance(nnode_, 0);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int flow = v_flow_[i];
      node_balance[s] += flow;
      node_balance[t] -= flow;
    }
    for (int i = 0; i < nnode_; ++i) {
      if (node_balance[i] > 0) {
        maxflow += node_balance[i];
      }
      // TODO: the imbalance needs to be accounted for
    }
    return maxflow;
  }

  long getMinCutValue() const { return mincut_value_; }

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

  void setMinCutSolution(const std::vector<bool> &new_solution) {
    std::copy(new_solution.begin(), new_solution.end(), x_.begin());
    computeMinCutValueInitial();
  }

  void setMinCutSolution(const std::vector<int> &new_solution) {
    std::copy(new_solution.begin(), new_solution.end(), x_.begin());
    computeMinCutValueInitial();
  }

private:
  void computeMinCutValueInitial() {
    mincut_value_ = 0;
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      if (x_[s] == 0 && x_[t] == 1) {
        mincut_value_ += forward_capacity;
      } else if (x_[s] == 1 && x_[t] == 0) {
        mincut_value_ += backward_capacity;
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      auto terminal_capacity = terminal_capacities_[i];
      if (x_[i] == 0 && terminal_capacity < 0) {
        mincut_value_ += -terminal_capacity;
      } else if (x_[i] == 1 && terminal_capacity > 0) {
        mincut_value_ += terminal_capacity;
      }
    }
    // add dual decomposition node potential terms (when/if applicable)
    size_t i = 0;
    for (const auto &constraint : dual_decomposition_constraints_) {
      int lagrange_multiplier_term = 0;
      for (const auto &arc_reference : constraint.source_arc_references) {
        lagrange_multiplier_term -= arc_reference->alpha;
      }
      for (const auto &arc_reference : constraint.target_arc_references) {
        lagrange_multiplier_term += arc_reference->alpha;
      }
      mincut_value_ +=
          lagrange_multiplier_term * x_[dual_decomposition_local_indices_[i++]];
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

  std::pair<int, int> arcGradients(int forward_capacity, int backward_capacity,
                                   int flow) const {
    int pos = flow + forward_capacity;
    int neg = flow - backward_capacity;
    return {pos, neg};
  }

  int nodeGradient(int terminal_capacity, int flow) const {
    return flow + terminal_capacity;
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
      int new_flow = 0;
      int dfp, dfn;
      if (pos < 0 || neg > 0) {
        dfp = std::min(pos, 0);
        dfn = std::max(neg, 0);
        if (std::abs(dfn) > std::abs(dfp)) {
          new_flow = -dfn;
        } else {
          new_flow = -dfp;
        }
      }
      if (new_flow != 0) {
        v_flow_[i] += new_flow;
        d_flow_[s] += new_flow;
        d_flow_[t] -= new_flow;
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

  void updateNodeTerminal(int i, int pos, bool do_update) {
    if (do_update) {
      auto existing_pos = maxflow_graph_.get_trcap(i);
      pos += existing_pos;
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
      auto pos = nodeGradient(terminal_capacities_[i], d_flow_[i]);
      updateNodeTerminal(i, pos, false);
    }
    size_t cache_index = 0;
    for (const auto &i : dual_decomposition_local_indices_) {
      auto pos = nodeGradient(terminal_capacities_[i], d_flow_[i]);
      pos += cached_lagrange_multipliers_[cache_index++];
      updateNodeTerminal(i, pos, false);
    }
  }

  void updateDualDecompositionNodePotentials() {
    // add dual decomposition node potential terms (when/if applicable)
    size_t cache_index = 0;
    for (const auto &index : dual_decomposition_local_indices_) {
      int pos = cached_lagrange_multipliers_[cache_index++];
      updateNodeTerminal(index, pos, true);
    }
  }

  void updateMinCut() {
    if (is_first_iteration_) {
      updateMinCutInitial();
    } else {
      updateMinCutIncremental();
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
      int lagrange_multiplier_term = cached_lagrange_multipliers_[cache_index];
      int last_lagrange_multiplier_term =
          cached_last_lagrange_multipliers_[cache_index++];
      if (last_lagrange_multiplier_term == lagrange_multiplier_term &&
          x_i_new == x_[index]) {
        continue;
      }
      mincut_value_ -= last_lagrange_multiplier_term * (x_[index]);
      mincut_value_ += lagrange_multiplier_term * (x_i_new);
    }

    // update node and arc terms that may have changed
    std::unordered_set<int> proccessed_arcs;
    auto nodes = maxflow_graph_.get_nodes();
    MaxflowGraph::arc_id first_arc = maxflow_graph_.get_first_arc();
    for (const int i : incremental_mincut_nodes_) {
      auto x_i_new =
          maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
      if (x_i_new == x_[i]) {
        continue; // do nothing
      }

      // process terminals
      auto terminal_capacity = terminal_capacities_[i];
      if (x_[i] == 0 && x_i_new == 1) {
        mincut_value_ += terminal_capacity;
      }
      if (x_[i] == 1 && x_i_new == 0) {
        mincut_value_ -= terminal_capacity;
      }

      // processes each possible arc
      MaxflowGraph::arc_id a;
      const auto &node_i = nodes[i];
      for (a = node_i.first; a; a = a->next) {
        auto arc_index = std::distance(first_arc, a) / 2;
        auto [iter, success] = proccessed_arcs.insert(arc_index);
        if (!success) {
          continue; // skip arc as it was processed with other node
        }
        auto forward_capacity = arc_capacities_[2 * arc_index + 0];
        auto backward_capacity = arc_capacities_[2 * arc_index + 1];
        int s = arcs_[2 * arc_index + 0];
        int t = arcs_[2 * arc_index + 1];
        auto x_s_new =
            maxflow_graph_.what_segment(s) == MaxflowGraph::SINK ? 1 : 0;
        auto x_t_new =
            maxflow_graph_.what_segment(t) == MaxflowGraph::SINK ? 1 : 0;
        if ((x_[s] == 0 && x_[t] == 1) && !(x_s_new == 0 && x_t_new == 1)) {
          mincut_value_ -= forward_capacity;
        }
        if ((x_[s] == 1 && x_[t] == 0) && !(x_s_new == 1 && x_t_new == 0)) {
          mincut_value_ -= backward_capacity;
        }
        if (!(x_[s] == 0 && x_[t] == 1) && (x_s_new == 0 && x_t_new == 1)) {
          mincut_value_ += forward_capacity;
        }
        if (!(x_[s] == 1 && x_[t] == 0) && (x_s_new == 1 && x_t_new == 0)) {
          mincut_value_ += backward_capacity;
        }
      }

      // update
      x_[i] = x_i_new;
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
      int new_flow = maxflow_graph_.get_rcap(a) - pos;
      v_flow_[i] += new_flow;
      d_flow_[s] += new_flow;
      d_flow_[t] -= new_flow;
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
      int new_flow = maxflow_graph_.get_rcap(first_arc + 2 * i) - pos;
      v_flow_[i] += new_flow;
      d_flow_[s] += new_flow;
      d_flow_[t] -= new_flow;
    }
  }

  void updateFlow() {
    if (is_first_iteration_) {
      updateFlowInitial();
    } else {
      updateFlowIncremental();
    }
  }

  void computeMaxflow() {
    if (is_first_iteration_) {
      maxflow_graph_.maxflow();
    } else {
      incremental_arcs_.clear();
      incremental_mincut_nodes_.clear();
      std::unordered_set<MaxflowGraph::arc_id> changed_arcs;
      maxflow_graph_.maxflow(true, changed_arcs, &maxflow_changed_list_);

      // update incremental nodes
      MaxflowGraph::node_id *ptr;
      for (ptr = maxflow_changed_list_.ScanFirst(); ptr;
           ptr = maxflow_changed_list_.ScanNext()) {
        MaxflowGraph::node_id i = *ptr;
        maxflow_graph_.remove_from_changed_list(i);
        incremental_mincut_nodes_.emplace_back(i);
      }
      maxflow_changed_list_.Reset();

      // update incremental arcs
      auto first_arc = maxflow_graph_.get_first_arc();
      for (const auto &arc_id : changed_arcs) {
        auto arc_index = std::distance(first_arc, arc_id) / 2;
        incremental_arcs_.emplace_back(arc_index);
      }
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
      int lagrange_multiplier_term = 0;
      int last_lagrange_multiplier_term = 0;
      for (const auto &arc_reference : constraint.source_arc_references) {
        lagrange_multiplier_term -= arc_reference->alpha;
        last_lagrange_multiplier_term -= arc_reference->last_alpha;
      }
      for (const auto &arc_reference : constraint.target_arc_references) {
        lagrange_multiplier_term += arc_reference->alpha;
        last_lagrange_multiplier_term += arc_reference->last_alpha;
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
  }

  /**
   * data passed into solver
   */
  int nnode_;
  int narc_;
  std::vector<int> arcs_;
  std::vector<int> arc_capacities_;
  std::vector<int> terminal_capacities_;

  /**
   * data structures needed for solving primal dual problem
   */
  std::vector<int> v_flow_; // flow on the arcs
  std::vector<int> d_flow_; // flow on the nodes
  std::vector<int> x_;      // mincut solution
  using MaxflowGraph =
      Graph</*captype=*/int, /*tcaptype=*/int, /*flowtype=*/long>;
  MaxflowGraph maxflow_graph_; // graph used to compute maxflow
  bool is_first_iteration_;
  bool is_first_iteration_of_new_scale_;

  Block<MaxflowGraph::node_id> maxflow_changed_list_;
  std::list<int> incremental_mincut_nodes_;
  std::list<int> incremental_arcs_;
  long mincut_value_;

  /**
   * specific to dual decomposition
   */
  struct DualDecompositionConstraint {
    std::vector<DualDecompositionConstraintArcReference> source_arc_references;
    std::vector<DualDecompositionConstraintArcReference> target_arc_references;
  };
  std::vector<int> dual_decomposition_local_indices_;
  std::vector<DualDecompositionConstraint> dual_decomposition_constraints_;

  std::vector<int> cached_lagrange_multipliers_;
  std::vector<int> cached_last_lagrange_multipliers_;
};

} // namespace mcpd3
