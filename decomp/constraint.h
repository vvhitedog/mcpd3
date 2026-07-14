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

#include <capacity.h>

#include <list>
#include <vector>

namespace mcpd3 {

struct DualDecompositionConstraintArc {
  Lagrange alpha; /* lagrange multiplier */
  Lagrange
      last_alpha; /* last lagrange multiplier recorded for incremental update */
  float alpha_momentum;       /* lagrange multiplier momentum */
  int partition_index_source; /* partition index for source node */
  int partition_index_target; /* partition index for target node */
  int local_index_source;     /* index within sub-problem of source */
  int local_index_target;     /* index within sub-problem of target */

  DualDecompositionConstraintArc(Lagrange alpha, Lagrange last_alpha,
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

struct DualDecompositionConstraintSnapshot {
  int constraint_id = -1;
  int global_node_id = -1;
  int partition_index_source = -1;
  int partition_index_target = -1;
  int local_index_source = -1;
  int local_index_target = -1;
  Lagrange alpha = 0;
  Lagrange last_alpha = 0;
  float alpha_momentum = 0;
};

struct DualDecompositionPartitionSnapshot {
  int partition_id = -1;
  Objective lower_bound = 0;
  Objective regularization_budget = 0;
  Objective regularization_contribution = 0;
  long regularization_anchor_sink_count = 0;
  long regularization_active_sink_count = 0;
  std::vector<int> local_labels;
};

using DualDecompositionConstraintArcReference =
    std::list<DualDecompositionConstraintArc>::iterator;

} // namespace mcpd3
