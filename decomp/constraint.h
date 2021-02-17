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

struct DualDecompositionConstraintArc {
  long alpha; /* lagrange multiplier */
  long
      last_alpha; /* last lagrange multiplier recorded for incremental update */
  float alpha_momentum;       /* lagrange multiplier momentum */
  int partition_index_source; /* partition index for source node */
  int partition_index_target; /* partition index for target node */
  int local_index_source;     /* index within sub-problem of source */
  int local_index_target;     /* index within sub-problem of target */

  DualDecompositionConstraintArc(int alpha, int last_alpha,
                                 float alpha_momentum,
                                 int partition_index_source,
                                 int partition_index_target,
                                 int local_index_source, int local_index_target)
      : alpha(alpha), last_alpha(last_alpha), alpha_momentum(alpha_momentum),
        partition_index_source(partition_index_source),
        partition_index_target(partition_index_target),
        local_index_source(local_index_source),
        local_index_target(local_index_target) {}
};

using DualDecompositionConstraintArcReference =
    std::list<DualDecompositionConstraintArc>::iterator;

} // namespace mcpd3
