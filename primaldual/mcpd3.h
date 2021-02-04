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
  PrimalDualMinCutSolver(int nnode, int narc, std::vector<int> &&arcs,
                         std::vector<int> arc_capacities,
                         std::vector<int> terminal_capacities)
      : nnode_(nnode), narc_(narc), arcs_(std::move(arcs)),
        arc_capacities_(std::move(arc_capacities)),
        terminal_capacities_(std::move(terminal_capacities)), v_flow_(narc_, 0),
        d_flow_(nnode_, 0), x_(nnode_, 0), maxflow_graph_(nnode_, narc_), is_first_iteration_(true), scale_(1), decomposition_size_(0) {
        //previous_terminal_node_capacities_(nnode_,-1) {
    initializeMaxflowGraph();
  }

  PrimalDualMinCutSolver(MinCutGraph min_cut_graph)
      : PrimalDualMinCutSolver(min_cut_graph.nnode, min_cut_graph.narc,
                               std::move(min_cut_graph.arcs),
                               std::move(min_cut_graph.arc_capacities),
                               std::move(min_cut_graph.terminal_capacities)) {}

  void setDecompositionSize(int decomposition_size) {
    decomposition_size_ = decomposition_size;
  }

  template<int scale>
  void scaleProblem() {
    scale_ *= scale;
    for (int i = 0; i < nnode_; ++i) {
      int &source_capacity = terminal_capacities_[2 * i + 0];
      int &sink_capacity = terminal_capacities_[2 * i + 1];
      source_capacity *= scale;
      sink_capacity *= scale;
      int &flow = d_flow_[i];
      flow *= scale;
    }
    for (int i = 0; i < narc_; ++i) {
      int &forward_capacity = arc_capacities_[2 * i + 0];
      int &backward_capacity = arc_capacities_[2 * i + 1];
      forward_capacity *= scale;
      backward_capacity *= scale;
      int &flow = v_flow_[i];
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
      maxflow_graph_.add_tweights(i, source_capacity, sink_capacity);
    }
    auto maxflow = maxflow_graph_.maxflow();
    maxflow_graph_.reset();
    initializeMaxflowGraph();
    return maxflow;
  }

  void solve() {
    //std::fill(x_.begin(),x_.end(),0); // reset min cut
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
    //if ( !is_first_iteration_ ) {
    //  markNodes();
    //}
    maxflow_graph_.maxflow(!is_first_iteration_);
    //maxflow_graph_.maxflow();
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
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      maxflow += std::min(source_capacity, sink_capacity);
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
      //auto [source_lagrange_multiplier,
      //target_lagrange_multiplier] = getDualDecompositionLagrangeMultiplier(i);
      //min_cut_value +=
      //    (source_lagrange_multiplier + target_lagrange_multiplier) * x_[i];
    }
    // add dual decomposition node potential terms (when/if applicable)
    for (const auto &[index,constraint] : dual_decomposition_constraints_map_ ) {
      int lagrange_multiplier_term = 0;
      if ( constraint.source_arc_reference ) {
        lagrange_multiplier_term -= (*constraint.source_arc_reference)->alpha;
      } 
      if ( constraint.target_arc_reference ) {
        lagrange_multiplier_term += (*constraint.target_arc_reference)->alpha;
      } 
      min_cut_value +=
          lagrange_multiplier_term * x_[index];
    }
    return min_cut_value;
  }

  //void addDualDecompositionConstraint(
  //    DualDecompositionConstraintArcReference arc_reference, bool is_source) {
  //  if (is_source) { // is source
  //    dual_decomposition_constraints_map_.insert(
  //        {arc_reference->local_index_source, {arc_reference, is_source}});
  //  } else { // is target
  //    dual_decomposition_constraints_map_.insert(
  //        {arc_reference->local_index_target, {arc_reference, is_source}});
  //  }
  //}

  void addSourceDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_source;
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if ( find_iter != dual_decomposition_constraints_map_.end() ) {
      auto &constraint = find_iter->second;
      constraint.source_arc_reference = arc_reference;
    } else {
      DualDecompositionConstraint constraint;
      constraint.source_arc_reference = arc_reference;
      dual_decomposition_constraints_map_.insert({index,constraint});
    }
  }

  void addTargetDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) {
    auto index = arc_reference->local_index_target;
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if ( find_iter != dual_decomposition_constraints_map_.end() ) {
      auto &constraint = find_iter->second;
      constraint.target_arc_reference = arc_reference;
    } else {
      DualDecompositionConstraint constraint;
      constraint.target_arc_reference = arc_reference;
      dual_decomposition_constraints_map_.insert({index,constraint});
    }
  }

  int getMinCutSolution( int index ) const {
    return x_[index];
  }

  void setMinCutSolution( const std::vector<int> new_solution ) {
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

  std::pair<int, int> arcGradients(int forward_capacity, int backward_capacity,
                                   int flow) const {
    int pos = flow + forward_capacity;
    int neg = flow - backward_capacity;
    return {pos, neg};
  }

  int nodeGradient(int source_capacity, int sink_capacity, int flow) const {
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

  void markNodes() {
    for (const auto &[index,constraint]: dual_decomposition_constraints_map_) {
      //if ( constraint.arc_reference->has_changed || 
      //    flow_modified_nodes_.find(index) != flow_modified_nodes_.end() ) {
      //if ( to_mark_nodes_.find(index) != to_mark_nodes_.end() ) {
        //maxflow_graph_.mark_node(index);
      //}
      //}
    }
  }

  std::pair<int,std::list<int>> getRegularizationParameters() const {
    std::list<int> to_regularize_indices;
    int lambda = 0;
    for (const auto &[index,constraint] : dual_decomposition_constraints_map_ ) {
      bool is_index_unsatisfied = false;
      if ( constraint.source_arc_reference ) {
        if ( (*constraint.source_arc_reference)->is_unsatisfied ){
          is_index_unsatisfied = true;
        }
      } 
      if ( constraint.target_arc_reference ) {
        if ( (*constraint.target_arc_reference)->is_unsatisfied ){
          is_index_unsatisfied = true;
        }
      } 
      if ( is_index_unsatisfied ) {
        to_regularize_indices.emplace_back(index);
        lambda += 1;
      }
    }
    lambda = lambda > 0 && decomposition_size_ > 0 ? scale_ / decomposition_size_ / lambda : 0;
    if ( lambda >= 1 ) {
      return {lambda,to_regularize_indices};
    }
    return {0,{}};
  }

  std::pair<int,int> getDualDecompositionLagrangeMultiplier(int index) const {
    std::pair<int,int> lagrange_multiplier_terms = {0,0};
    auto constraint_map_iter = dual_decomposition_constraints_map_.find(index);
    if (constraint_map_iter != dual_decomposition_constraints_map_.end()) {
      const auto &constraint = constraint_map_iter->second;
      if ( constraint.source_arc_reference ) {
        lagrange_multiplier_terms.first -= (*constraint.source_arc_reference)->alpha;
      }
      if ( constraint.target_arc_reference ) {
        lagrange_multiplier_terms.first += (*constraint.target_arc_reference)->alpha;
      }
      //if (constraint.is_source) { // is source in constraint
      //  lagrange_multiplier_term -= constraint.arc_reference->alpha;
      //} else { // is target in constraint
      //  lagrange_multiplier_term += constraint.arc_reference->alpha;
      //}
    }
    return lagrange_multiplier_terms;
  }
  
  void updateNodeTerminal(int i, int pos, int &nviolated, bool do_update) {
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
      int source_capacity = terminal_capacities_[2 * i + 0];
      int sink_capacity = terminal_capacities_[2 * i + 1];
      auto pos = nodeGradient(source_capacity, sink_capacity, d_flow_[i]);
      updateNodeTerminal(i,pos,nviolated,false);
      // include any active dual decomposition constraint if applicable
      //auto [source_lagrange_multiplier,
      //target_lagrange_multiplier] = getDualDecompositionLagrangeMultiplier(i);
      //pos += source_lagrange_multiplier + target_lagrange_multiplier;
      //if (pos < 0) {
      //  nviolated++;
      //}
      //if ( !is_first_iteration_ ) {
      //  auto stored_pos = maxflow_graph_.get_trcap(i);
      //  if ( stored_pos != pos ) {
      //    maxflow_graph_.set_trcap(i, pos);
      //    maxflow_graph_.mark_node(i);
      //  }
      //} else if (is_first_iteration_) {
      //  maxflow_graph_.set_trcap(i, pos);
      //}
    }
    // add dual decomposition node potential terms (when/if applicable)
    for (const auto &[index,constraint] : dual_decomposition_constraints_map_ ) {
      int pos = 0;
      if ( constraint.source_arc_reference ) {
        pos -= (*constraint.source_arc_reference)->alpha;
      } 
      if ( constraint.target_arc_reference ) {
        pos += (*constraint.target_arc_reference)->alpha;
      } 
      updateNodeTerminal(index,pos,nviolated,true);
    }
    // add regularization node potential terms (when/if applicable)
    auto [lambda,to_regularize_indices] = getRegularizationParameters();
    for(auto i : to_regularize_indices ) {
      updateNodeTerminal(i,lambda,nviolated,true);
      //auto pos = maxflow_graph_.get_trcap(i);
      //pos += lambda;
      //if ( pos < 0 ) {
      //  nviolated++;
      //}
      //if ( !is_first_iteration_ ) {
      //  maxflow_graph_.set_trcap(i, pos);
      //  maxflow_graph_.mark_node(i);
      //} else if (is_first_iteration_) {
      //  maxflow_graph_.set_trcap(i, pos);
      //}
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
      auto [pos, neg] = arcGradients(forward_capacity, backward_capacity, flow);
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
  bool is_first_iteration_;
  int scale_;
  int decomposition_size_;
  //std::unordered_set<int> to_mark_nodes_;
  //std::vector<int> previous_terminal_node_capacities_; // from a previous iteration

  /**
   * specific to dual decomposition
   */
  struct DualDecompositionConstraint {
    std::optional<DualDecompositionConstraintArcReference> source_arc_reference;
    std::optional<DualDecompositionConstraintArcReference> target_arc_reference;
    DualDecompositionConstraint(): source_arc_reference({}),
      target_arc_reference({}) {}
  };
  std::unordered_map</*local_index=*/int, DualDecompositionConstraint>
      dual_decomposition_constraints_map_;
};

} // namespace mcpd3
