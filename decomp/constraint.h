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

#include <list>

namespace mcpd3 {

  constexpr int MAX_DUAL_DECOMPOSITION_RECURSION_LEVEL = 10;

struct DualDecompositionConstraintArc {
  long alpha; /* lagrange multiplier */
  long
      last_alpha; /* last lagrange multiplier recorded for incremental update */
  float alpha_momentum;       /* lagrange multiplier momentum */
  int partition_index_source[MAX_DUAL_DECOMPOSITION_RECURSION_LEVEL]; /* partition index for source node */
  int partition_index_target[MAX_DUAL_DECOMPOSITION_RECURSION_LEVEL]; /* partition index for target node */
  int local_index_source[MAX_DUAL_DECOMPOSITION_RECURSION_LEVEL];     /* index within sub-problem of source */
  int local_index_target[MAX_DUAL_DECOMPOSITION_RECURSION_LEVEL];     /* index within sub-problem of target */

  DualDecompositionConstraintArc(int level, 
                                 int _partition_index_source,
                                 int _partition_index_target,
                                 int _local_index_source,
                                 int _local_index_target)
      : alpha(0), last_alpha(0), alpha_momentum(0)
        {
          partition_index_source[level] = _partition_index_source;
          partition_index_target[level] = _partition_index_target;
          local_index_source[level] = _local_index_source;
          local_index_target[level] = _local_index_target;
        }
};

using DualDecompositionConstraintArcReference =
    std::list<DualDecompositionConstraintArc>::iterator;

} // namespace mcpd3
