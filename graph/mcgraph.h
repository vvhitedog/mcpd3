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

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcpd3 {

struct MinCutGraph {
  int nnode, narc;
  std::vector<int> arcs;
  std::vector<int> arc_capacities;
  std::vector<int> terminal_capacities;
};

void to_file(const std::string &filename, const MinCutGraph &min_cut_graph) {
  FILE *stream = fopen(filename.c_str(), "w");
  if (!stream) {
    throw std::runtime_error("Failed to open filename '" + filename +
                             "' for writing.");
  }
  fprintf(stream, "n %d\n", min_cut_graph.nnode);
  fprintf(stream, "m %d\n", min_cut_graph.narc);
  for (int i = 0; i < min_cut_graph.narc; ++i) {
    int s = min_cut_graph.arcs[2 * i + 0];
    int t = min_cut_graph.arcs[2 * i + 1];
    int forward_capacity = min_cut_graph.arc_capacities[2 * i + 0];
    int backward_capacity = min_cut_graph.arc_capacities[2 * i + 1];
    fprintf(stream, "a %8d %8d %8d %8d\n", s, t, forward_capacity,
            backward_capacity);
  }
  for (int i = 0; i < min_cut_graph.nnode; ++i) {
    int source_capacity = min_cut_graph.terminal_capacities[2 * i + 0];
    int sink_capacity = min_cut_graph.terminal_capacities[2 * i + 1];
    fprintf(stream, "t %8d %8d %8d\n", i, source_capacity, sink_capacity);
  }
  fclose(stream);
}

} // namespace mcpd3
