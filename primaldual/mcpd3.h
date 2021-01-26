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

#include <iostream>

namespace mcpd3 {

// TODO: erase me
// ambig - forward
// cost - backward
// source - ambig
// sink - cost

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

  void solve() {
    resetSolver();
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
    auto mfval = maxflow_graph_.maxflow();
    std::cout << "computed internal mf val: " << mfval << "\n";
    updateFlow();   // get updated flow
    updateMinCut(); // get updated min cut solution
  }

  std::vector<int>& terminal_capacities() { 
    return terminal_capacities_;
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

private:
  void resetSolver() {
    //static int c = 0;
    //std::srand(time(NULL) + c++);
    //for (int i = 0; i < nnode_; ++i) {
    //  x_[i] = std::rand() % 2;
    //  //x_[i] = 1;
    //  assert(x_[i] == 0 || x_[i] == 1);
    //}

    //std::fill(x_.begin(), x_.end(), 0);
    //x_ = {1,0};

    // XXX: interesting to comment out below to check to see if
    // complementary slackness conditions can be used to find maxflow from
    // mincut
    //
    //std::fill(v_flow_.begin(), v_flow_.end(), 0);
    //std::fill(d_flow_.begin(), d_flow_.end(), 0);

     for (int i = 0; i < nnode_; ++i) {
      // NOTE: resetting terminal capacities to 0 is necessary as
      // `add_tweights` is used below
      maxflow_graph_.set_trcap(i, 0);
    }

    maxflow_graph_.reset();
    initializeMaxflowGraph();

  }

  void initializeMaxflowGraph() {
    maxflow_graph_.add_node(nnode_);
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      maxflow_graph_.add_edge(s, t, 0, 0);
    }
  }

  std::pair<int, int> arcGradients(int x_s, int x_t, int forward_capacity,
                                   int backward_capacity, int flow) {
    int pos = flow;
    int neg = flow;
    pos += forward_capacity *
           (std::max(x_t - x_s + 1, 0) - std::max(x_t - x_s, 0));
    neg += forward_capacity *
           (std::max(x_t - x_s, 0) - std::max(x_t - x_s - 1, 0));
    pos += backward_capacity *
           (std::max(x_s - x_t - 1, 0) - std::max(x_s - x_t, 0));
    neg += backward_capacity *
           (std::max(x_s - x_t, 0) - std::max(x_s - x_t + 1, 0));
    return {pos, neg};
  }

  std::pair<int,int> nodeGradient(int x_s, int source_capacity, int sink_capacity, int flow) {
    // NOTE: the negative gradient can be calculated as:
    //
    // int neg = flow;
    // neg += source_capacity * (std::max(x_s,0) - std::max(x_s - 1,0));
    // neg += sink_capacity * (std::max(1 - x_s,0) - std::max(1 - x_s + 1,0));
    //
    // but there is no point, it is not used.
    int pos = flow;
    pos += source_capacity * (std::max(x_s + 1, 0) - std::max(x_s, 0));
    pos += sink_capacity * (std::max(1 - x_s - 1, 0) - std::max(1 - x_s, 0));
    int neg = flow;
    neg += sink_capacity * (std::max(x_s,0) - std::max(x_s - 1,0));
    neg += source_capacity * ( std::max(1 - x_s,0) - std::max(1 - x_s + 1,0));
    //neg += sink_capacity * (std::max(x_s,0) - std::max(x_s - 1,0));
    //neg += source_capacity * (std::max(1 - x_s,0) - std::max(1 - x_s - 1,0));
    return {pos,neg};
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
          arcGradients(x_[s], x_[t], forward_capacity, backward_capacity, flow);
          //arcGradients(0, 0, forward_capacity, backward_capacity, flow);
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
      std::tie(pos, neg) = arcGradients(x_[s], x_[t], forward_capacity,
                                        backward_capacity, v_flow_[i]);
      //std::tie(pos, neg) = arcGradients(0, 0, forward_capacity,
      //                                  backward_capacity, v_flow_[i]);
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
      auto [pos,neg] =
          nodeGradient(x_[i], source_capacity, sink_capacity, d_flow_[i]);
          //nodeGradient(0, source_capacity, sink_capacity, d_flow_[i]);

      if ( x_[i] == 0 ) {
        if ( pos < 0 ) {
          nviolated++;
        }
        std::cout << "attaching pos: " << pos << " to node " << i << std::endl;
        maxflow_graph_.set_trcap(i, pos);
      } else if ( x_[i] == 1 ) {
        if ( neg > 0 ) {
          nviolated++;
        }
        maxflow_graph_.set_trcap(i, neg);
        std::cout << "attaching neg: " << neg << " to node " << i << std::endl;
      }

      // if (pos > 0) {
      //  maxflow_graph_.add_tweights(i, pos, 0);
      //} else {
      //  maxflow_graph_.add_tweights(i, 0, -pos);
      //}
    }
    return nviolated;
  }

  void updateMinCut() {
    for (int i = 0; i < nnode_; ++i) {
      std::cout << " decoded : " << maxflow_graph_.what_segment(i) << " for node " << i << " which had init soln: " << x_[i] << "\n";
      if ( x_[i] == 0 ) {
      x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 1 : 0;
      } else if (x_[i] == 1) {
      x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK ? 0 : 1;
      }
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
          arcGradients(x_[s], x_[t], forward_capacity, backward_capacity, flow);
          //arcGradients(0, 0, forward_capacity, backward_capacity, flow);
      //int32_t new_flow = maxflow_graph_.get_rcap(a) - pos;
      int new_flow = maxflow_graph_.get_rcap(a) - pos;
      v_flow_[i] += new_flow;
      d_flow_[s] += new_flow;
      d_flow_[t] -= new_flow;
      a = maxflow_graph_.get_next_arc(a);
      a = maxflow_graph_.get_next_arc(a);
    }
  }

  /**
   * Data passed into solver.
   */
  int nnode_;
  int narc_;
  std::vector<int> arcs_;
  std::vector<int> arc_capacities_;
  std::vector<int> terminal_capacities_;

  /**
   * Data structures needed for solving primal dual problem.
   */
  std::vector<int> v_flow_; // flow on the arcs
  std::vector<int> d_flow_; // flow on the nodes
  std::vector<int> x_;      // mincut solution
  using MaxflowGraph =
      Graph</*captype=*/int, /*tcaptype=*/int, /*flowtype=*/long>;
  MaxflowGraph maxflow_graph_; // graph used to compute maxflow
};

} // namespace mcpd3
