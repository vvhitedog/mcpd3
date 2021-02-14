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
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>

#include <decomp/constraint.h>
#include <graph/mcgraph.h>
#include <maxflow/graph.h>

namespace mcpd3 {

class PrimalDualMinCutSolver {
public:
  PrimalDualMinCutSolver(int nnode, int narc, std::vector<long> &&arcs,
                         std::vector<long> arc_capacities,
                         std::vector<long> terminal_capacities)
      : nnode_(nnode), narc_(narc), arcs_(std::move(arcs)),
        arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)), v_flow_(narc_, 0),
        d_flow_(nnode_, 0), x_(nnode_, 0), maxflow_graph_(nnode_, narc_), is_first_iteration_(true) {
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

    std::vector<bool> visited(nnode_,false);
    std::vector<int> dist(nnode_,0);
    std::list<int> q;

    int index = 0;
    for ( const auto &index : seeds ) {
      q.emplace_back(index);
      visited[index] = true;
    }
    MaxflowGraph::arc_id a;
    auto nodes = maxflow_graph_.get_nodes();


    while( !q.empty() ) {
      int u = q.front();
      q.pop_front();
      if (dist[u] >= rad) {
        continue;
      }
      const auto &_u = nodes[u];
			for (a=_u.first; a; a=a->next) {
        auto v = a->head;
        auto iv = std::distance(nodes,v);
        if ( !visited[iv] ) {
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
      if ( !visited[s] || !visited[t] ) {
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
      if ( !visited[i] ) {
        maxflow_graph_.set_trcap(i, 0);
        continue;
      }
      num_nodes_in_decoding++;
      if ( dist[i] < rad ) {
        auto source_capacity = terminal_capacities_[2 * i + 0];
        auto sink_capacity = terminal_capacities_[2 * i + 1];
        maxflow_graph_.add_tweights(i, source_capacity, sink_capacity);
      } else {
        if ( x_[i] == 0 ) {
          maxflow_graph_.set_trcap(i, std::numeric_limits<int>::max()/2);
        } else {
          maxflow_graph_.set_trcap(i, -std::numeric_limits<int>::max()/2);
        }
      }
    }

    auto maxflow_val = maxflow_graph_.maxflow();

    for (int i = 0; i < nnode_; ++i) {
      if ( !visited[i] ) {
        continue;
      }
      if ( dist[i] < rad ) {
        x_[i] = maxflow_graph_.what_segment(i) == MaxflowGraph::SINK;
      }
    }

    printf("/////////////////////////////////////////\n");
    printf("//////////  DECODING STATS //////////////\n");
    printf("// num_arcs in decoding :  %8lu    //\n",num_arcs_in_decoding);
    printf("// num_nodes in decoding : %8lu    //\n",num_nodes_in_decoding);
    printf("// maxflow val :           %8ld    //\n",maxflow_val);
    printf("/////////////////////////////////////////\n");

  }


  template<int scale>
  void scaleProblem() {
    for (int i = 0; i < nnode_; ++i) {
      auto &source_capacity = terminal_capacities_[2 * i + 0];
      auto &sink_capacity = terminal_capacities_[2 * i + 1];
      source_capacity *= scale;
      sink_capacity *= scale;
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
      maxflow_graph_.set_rcap(a, flow*scale);
      a = maxflow_graph_.get_next_arc(a);
      flow = maxflow_graph_.get_rcap(a);
      maxflow_graph_.set_rcap(a, flow*scale);
      a = maxflow_graph_.get_next_arc(a);
    }
    for (int i = 0; i < nnode_; ++i) {
      int flow;
      flow = maxflow_graph_.get_trcap(i);
      maxflow_graph_.set_trcap(i,flow*scale);
    }
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
      auto source_capacity = terminal_capacities_[2 * i + 0];
      auto sink_capacity = terminal_capacities_[2 * i + 1];
      maxflow_graph_.add_tweights(i, source_capacity, sink_capacity);
    }
    auto maxflow = maxflow_graph_.maxflow();
    maxflow_graph_.reset();
    initializeMaxflowGraph();
    return maxflow;
  }

  void solve() {
    if ( is_first_iteration_ ) {
    initializeFlow(); // finds a flow satisfying arc based lagrange multiplier
                      // complementary slackness conditions
    }
    auto nviolated =
        updateNodePotentials(); // finds which node based lagrange multiplier
                                // complementary slackness conditions are
                                // violated and sets source and sink capacities
                                // accordingly
    if (!nviolated) {
      std::fill(x_.begin(),x_.end(),0); // reset min cut; no flow adjustment is necessary
      return; // no violated node based lagrange multiplier complementary
              // slackness conditions; translation: solution is optimal already
    }
    maxflow_graph_.maxflow(!is_first_iteration_);
    is_first_iteration_ = false;
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
    for (int i = 0; i < nnode_; ++i) {
      if (node_balance[i] > 0) {
        maxflow += node_balance[i];
      }
      auto source_capacity = terminal_capacities_[2 * i + 0];
      auto sink_capacity = terminal_capacities_[2 * i + 1];
      maxflow += std::min(source_capacity, sink_capacity);
    }
    return maxflow;
  }

  long getMinCutValue() const {
    long min_cut_value = 0;
    for (int i = 0; i < narc_; ++i) {
      int s = arcs_[2 * i + 0];
      int t = arcs_[2 * i + 1];
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      if (x_[s] == 0 && x_[t] == 1) {
        min_cut_value += forward_capacity;
      } else if (x_[s] == 1 && x_[t] == 0) {
        min_cut_value += backward_capacity;
      }
    }
    for (int i = 0; i < nnode_; ++i) {
      auto source_capacity = terminal_capacities_[2 * i + 0];
      auto sink_capacity = terminal_capacities_[2 * i + 1];
      if (x_[i] == 0) {
        min_cut_value += sink_capacity;
      } else {
        min_cut_value += source_capacity;
      }
    }
    // add dual decomposition node potential terms (when/if applicable)
    for (const auto &[index,constraint] : dual_decomposition_constraints_map_ ) {
      int lagrange_multiplier_term = 0;
      for (const auto &arc_reference : constraint.source_arc_references ) {
        lagrange_multiplier_term -= arc_reference->alpha;
      }
      for (const auto &arc_reference : constraint.target_arc_references ) {
        lagrange_multiplier_term += arc_reference->alpha;
      }
      min_cut_value +=
          lagrange_multiplier_term * x_[index];
    }
    return min_cut_value;
  }

  void addSourceDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_source;
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if ( find_iter != dual_decomposition_constraints_map_.end() ) {
      auto &constraint = find_iter->second;
      constraint.source_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.source_arc_references.emplace_back(arc_reference);
      dual_decomposition_constraints_map_.insert({index,constraint});
    }
  }

  void addTargetDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_target;
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if ( find_iter != dual_decomposition_constraints_map_.end() ) {
      auto &constraint = find_iter->second;
      constraint.target_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.target_arc_references.emplace_back(arc_reference);
      dual_decomposition_constraints_map_.insert({index,constraint});
    }
  }

  int getMinCutSolution( int index ) const {
    return x_[index];
  }

  void setMinCutSolution( const std::vector<long> new_solution ) {
    std::copy(new_solution.begin(),new_solution.end(),x_.begin());
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

  std::pair<long, long> arcGradients(long forward_capacity, long backward_capacity,
                                   long flow) const {
    long pos = flow + forward_capacity;
    long neg = flow - backward_capacity;
    return {pos, neg};
  }

  long nodeGradient(long source_capacity, long sink_capacity, long flow) const {
    return flow + source_capacity - sink_capacity;
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
      long new_flow = 0;
       long dfp, dfn;
      if (pos < 0 || neg > 0) {
        dfp = std::min(pos, 0L);
        dfn = std::max(neg, 0L);
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

  std::pair<long,long> getDualDecompositionLagrangeMultiplier(int index) const {
    std::pair<long,long> lagrange_multiplier_terms = {0,0};
    auto constraint_map_iter = dual_decomposition_constraints_map_.find(index);
    if (constraint_map_iter != dual_decomposition_constraints_map_.end()) {
      const auto &constraint = constraint_map_iter->second;
      for ( const auto &arc_reference : constraint.source_arc_references ) {
        lagrange_multiplier_terms.first -= arc_reference->alpha;
      }
      for ( const auto &arc_reference : constraint.target_arc_references ) {
        lagrange_multiplier_terms.first += arc_reference->alpha;
      }
    }
    return lagrange_multiplier_terms;
  }
  
  void updateNodeTerminal(int i, long pos, int &nviolated, bool do_update) {
    if (do_update ) {
        auto existing_pos = maxflow_graph_.get_trcap(i);
        if ( existing_pos < 0 ) {
          nviolated--;
        }
        pos += existing_pos;
        if ( pos < 0 ) {
          nviolated++;
        }
        if ( !is_first_iteration_ ) {
          if ( existing_pos != pos ) {
            maxflow_graph_.set_trcap(i, pos);
            maxflow_graph_.mark_node(i);
          }
        } else if (is_first_iteration_) {
          maxflow_graph_.set_trcap(i, pos);
        }
    } else {
      if (pos < 0) {
        nviolated++;
      }
      if ( !is_first_iteration_ ) {
        auto stored_pos = maxflow_graph_.get_trcap(i);
        if ( stored_pos != pos ) {
          maxflow_graph_.set_trcap(i, pos);
          maxflow_graph_.mark_node(i);
        }
      } else if (is_first_iteration_) {
        maxflow_graph_.set_trcap(i, pos);
      }
    }
  }

  int updateNodePotentials() {
    int nviolated = 0;
    // add mincut node potential terms
    for (int i = 0; i < nnode_; ++i) {
      auto source_capacity = terminal_capacities_[2 * i + 0];
      auto sink_capacity = terminal_capacities_[2 * i + 1];
      auto pos = nodeGradient(source_capacity, sink_capacity, d_flow_[i]);
      updateNodeTerminal(i,pos,nviolated,false);
    }
    // add dual decomposition node potential terms (when/if applicable)
    for (const auto &[index,constraint] : dual_decomposition_constraints_map_ ) {
      long pos = 0;
      for (const auto &arc_reference : constraint.source_arc_references ) {
        pos -= arc_reference->alpha;
      }
      for (const auto &arc_reference : constraint.target_arc_references ) {
        pos += arc_reference->alpha;
      }
      updateNodeTerminal(index,pos,nviolated,true);
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
      auto forward_capacity = arc_capacities_[2 * i + 0];
      auto backward_capacity = arc_capacities_[2 * i + 1];
      auto flow = v_flow_[i];
      auto [pos, neg] = arcGradients(forward_capacity, backward_capacity, flow);
      long new_flow = maxflow_graph_.get_rcap(a) - pos;
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
  std::vector<long> arcs_;
  std::vector<long> arc_capacities_;
  std::vector<long> terminal_capacities_;

  /**
   * data structures needed for solving primal dual problem
   */
  std::vector<long> v_flow_; // flow on the arcs
  std::vector<long> d_flow_; // flow on the nodes
  std::vector<long> x_;      // mincut solution
  using MaxflowGraph =
      Graph</*captype=*/long, /*tcaptype=*/long, /*flowtype=*/long>;
  MaxflowGraph maxflow_graph_; // graph used to compute maxflow
  bool is_first_iteration_;

  /**
   * specific to dual decomposition
   */
  struct DualDecompositionConstraint {
    std::list<DualDecompositionConstraintArcReference> source_arc_references;
    std::list<DualDecompositionConstraintArcReference> target_arc_references;
  };
  std::unordered_map</*local_index=*/int, DualDecompositionConstraint>
      dual_decomposition_constraints_map_;

};

} // namespace mcpd3
