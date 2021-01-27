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

#include <tuple>
#include <vector>

#include <maxflow/graph.h>
#include <graph/mcgraph.h>
#include <decomp/constraint.h>

namespace mcpd3 {

class PrimalDualMinCutSolver {
public:
  PrimalDualMinCutSolver(int nnode, int narc, std::vector<int> &&arcs,
                         std::vector<int> &&arc_capacities,
                         std::vector<int> &&terminal_capacities)
      : nnode_(nnode), narc_(narc), arcs_(arcs),
        arc_capacities_(arc_capacities),
        terminal_capacities_(terminal_capacities), v_flow_(narc_, 0),
        d_flow_(nnode_, 0), x_(nnode_, 0), maxflow_graph_(nnode_, narc_) {
    initializeMaxflowGraph();
  }

  PrimalDualMinCutSolver(MinCutGraph &min_cut_graph)
      : PrimalDualMinCutSolver(min_cut_graph.nnode, min_cut_graph.narc,
          std::move(min_cut_graph.arcs),
          std::move(min_cut_graph.arc_capacities), std::move(min_cut_graph.terminal_capacities) )  {}

  long maxflow() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      maxflow_graph_.set_rcap(a, forward_capacity);
      a = maxflow_graph_.get_next_arc(a);
      maxflow_graph_.set_rcap(a, backward_capacity);
      a = maxflow_graph_.get_next_arc(a);
    }
    for (int i = 0; i < nnode_; ++i) {
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      maxflow_graph_.add_tweights(i,source_capacity,sink_capacity);
    }
    auto maxflow = maxflow_graph_.maxflow();
    maxflow_graph_.reset();
    initializeMaxflowGraph();
    return maxflow;
  }

  void solve() {
    initializeFlow(); // finds a flow satisfying arc based lagrange multiplier
                      // complementary slackness conditions
    auto nviolated =
        updateNodePotentials(); // finds which node based lagrange multiplier
                                // complementary slackness conditions are
                                // violated and sets source and sink capacities
                                // accordingly
    if (!nviolated) {
      return; // no violated node based lagrange multiplier complementary
              // slackness conditions; translation: solution is optimal already
    }
    maxflow_graph_.maxflow();
    updateFlow();   // get updated flow
    updateMinCut(); // get updated min cut solution
  }

  long getMaxFlowValue() const {
    long maxflow = 0;
    std::vector<long> node_balance(nnode_, 0);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int flow = v_flow_[i];
      node_balance[s] += flow;
      node_balance[t] -= flow;
    }
    for ( int i = 0; i < nnode_; ++i ) {
      if ( node_balance[i] > 0 ) {
        maxflow += node_balance[i];
      }
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      maxflow += std::min(source_capacity,sink_capacity);
    }
    return maxflow;
  }

  long getMinCutValue() const {
    long min_cut_value = 0;
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      if (x_[s] == 0 && x_[t] == 1) {
        min_cut_value += forward_capacity;
      } else if (x_[s] == 1 && x_[t] == 0) {
        min_cut_value += backward_capacity;
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      if (x_[i] == 0) {
        min_cut_value += sink_capacity;
      } else {
        min_cut_value += source_capacity;
      }
    }
    return min_cut_value;
  }

  void addDualDecompositionConstraint( DualDecompositionConstraintArcReference arc_reference,
      bool is_source ) {
    dual_decomposition_constraints_.emplace_back(is_source,arc_reference);
  }

private:
  void initializeMaxflowGraph() {
    maxflow_graph_.add_node(nnode_);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      maxflow_graph_.add_edge(s, t, 0, 0);
    }
  }

  std::pair<int, int> arcGradients(int forward_capacity,
                                   int backward_capacity, int flow) const {
    int pos = flow + forward_capacity;
    int neg = flow - backward_capacity;
    return {pos, neg};
  }

  int nodeGradient( int source_capacity, int sink_capacity, int flow) const {
    return flow + source_capacity - sink_capacity;
  }

  void initializeFlow() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      int flow = v_flow_[i];
      auto [pos, neg] =
          arcGradients(forward_capacity, backward_capacity, flow);
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
      std::tie(pos, neg) = arcGradients(forward_capacity,
                                        backward_capacity, v_flow_[i]);
      assert(pos >= 0 && neg <= 0);
      maxflow_graph_.set_rcap(a, pos);
      a = maxflow_graph_.get_next_arc(a);
      maxflow_graph_.set_rcap(a, -neg);
      a = maxflow_graph_.get_next_arc(a);
    }
  }

  int updateNodePotentials() {
    int nviolated = 0;
    for (int i = 0; i < nnode_; ++i) {
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      auto pos =
          nodeGradient(source_capacity, sink_capacity, d_flow_[i]);
      if ( pos < 0 ) {
        nviolated++;
      }
      maxflow_graph_.set_trcap(i, pos);
    }
    return nviolated;
  }

  void updateMinCut() {
    for (int i = 0; i < nnode_; ++i) {
      x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
    }
  }

  void updateFlow() {
    MaxflowGraph::arc_id a = maxflow_graph_.get_first_arc();
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      int forward_capacity = arc_capacities_[2 * i + 0];
      int backward_capacity = arc_capacities_[2 * i + 1];
      int flow = v_flow_[i];
      auto [pos, neg] =
          arcGradients(forward_capacity, backward_capacity, flow);
      int new_flow = maxflow_graph_.get_rcap(a) - pos;
      v_flow_[i] += new_flow;
      d_flow_[s] += new_flow;
      d_flow_[t] -= new_flow;
      a = maxflow_graph_.get_next_arc(a);
      a = maxflow_graph_.get_next_arc(a);
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

  /**
   * specific to dual decomposition
   */
  struct DualDecompositionConstraint {
    bool is_source; // otherwise is_target
    DualDecompositionConstraintArcReference dual_decomposition_constraint_arc_ref;
  };
  std::list<DualDecompositionConstraint> dual_decomposition_constraints_;
};

} // namespace mcpd3
