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

#include <decomp/constraint.h>
#include <unordered_map>

namespace mcpd3 {

class DualDecompositionSubSolver {
public:
  ~DualDecompositionSubSolver() = default;

  virtual void addSourceDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) = 0;

  virtual void addTargetDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference) = 0;

  virtual void solve() = 0;

  virtual long getMinCutValue() const = 0;

  virtual int getMinCutSolution(int index) const = 0;

protected:
  void addSourceDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference, int index) {
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if (find_iter != dual_decomposition_constraints_map_.end()) {
      auto &constraint = find_iter->second;
      constraint.source_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.source_arc_references.emplace_back(arc_reference);
      dual_decomposition_constraints_map_.insert({index, constraint});
    }
  }

  void addTargetDualDecompositionConstraint(
      DualDecompositionConstraintArcReference arc_reference, int index) {
    auto find_iter = dual_decomposition_constraints_map_.find(index);
    if (find_iter != dual_decomposition_constraints_map_.end()) {
      auto &constraint = find_iter->second;
      constraint.target_arc_references.emplace_back(arc_reference);
    } else {
      DualDecompositionConstraint constraint;
      constraint.target_arc_references.emplace_back(arc_reference);
      dual_decomposition_constraints_map_.insert({index, constraint});
    }
  }

  /**
   * storage of dual decomposition constriants
   */
  struct DualDecompositionConstraint {
    std::list<DualDecompositionConstraintArcReference> source_arc_references;
    std::list<DualDecompositionConstraintArcReference> target_arc_references;
  };
  std::unordered_map</*local_index=*/int, DualDecompositionConstraint>
      dual_decomposition_constraints_map_;
};

} // namespace mcpd3
