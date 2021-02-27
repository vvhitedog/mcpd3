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

#include <cassert>
#include <string>
#include <unordered_map>

#include <graph/mcgraph.h>

namespace mcpd3 {

namespace _dimacs_implementation {

int remap_index(int original_index, int source_index, int sink_index) {
  int new_index = original_index;
  if (original_index > source_index) {
    --new_index;
  }
  if (original_index > sink_index) {
    --new_index;
  }
  return --new_index; // indices in DIMACS start at 1 but we wish for them
                      // to start at 0
}

template <typename ArcOperator, typename TerminalOperator>
void read_dimacs_general(const std::string &filename, 
                         ArcOperator arc_op, TerminalOperator term_op) {
  const int line_length = 1024;
  char line[line_length];
  int n, m, s, t, cap, source, sink, _s, _t;
  char c;
  FILE *stream = nullptr;
  source = -1;
  sink = -1;
  if (!(stream = fopen(filename.c_str(), "r"))) {
    std::string err_msg = "failed to open file for reading: " + filename;
    throw std::runtime_error(err_msg.c_str());
  }
  while (fgets(line, line_length, stream)) {
    switch (*line) {
    case 'p':
      if (sscanf(line, "%*c %*s %d %d", &n, &m) != 2) {
        std::string err_msg =
            "p line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      break;
    case 'a':
      if (source == -1 || sink == -1) {
        std::string err_msg =
            "'a' line occured beforce setting source/sink in DIMACS file:" +
            filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      if (sscanf(line, "%*c %d %d %d", &s, &t, &cap) != 3) {
        std::string err_msg =
            "'a' line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      if (t == source || s == sink) {
        std::string err_msg =
            "specified source or sink as target or source node incorrectly:" +
            filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      if (s != source && t != sink) { // handle non terminal arc
        _s = remap_index(s, source, sink);
        _t = remap_index(t, source, sink);
        arc_op(_s, _t, cap);
      } else { // handle terminal arc
        if (s == source) {
          _t = remap_index(t, source, sink);
          term_op(true, _t, cap);
        } else if (t == sink) {
          _s = remap_index(s, source, sink);
          term_op(false, _s, cap);
        }
      }
      break;
    case 'n':
      if (sscanf(line, "%*c %d %c", &s, &c) != 2) {
        std::string err_msg =
            "'n' line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      if (c == 's') {
        source = s;
      }
      if (c == 't') {
        sink = s;
      }
      break;
    case 'c':
      // comment
      break;
    default:
      break;
    }
  }
  fclose(stream);
}
} // namespace _dimacs_implementation

MinCutGraph read_dimacs(const std::string &filename) {
  MinCutGraph g;
  g.nnode = 0;
  std::unordered_map<int, std::unordered_map<int, int>> arc_adjacency;

  auto arc_op = [&](int s, int t, int cap) {
    g.nnode = std::max(g.nnode,s+1);
    g.nnode = std::max(g.nnode,t+1);
    if (arc_adjacency[s].find(t) == arc_adjacency[s].end()) {
      arc_adjacency[s][t] = cap;
      if (arc_adjacency[t].find(s) == arc_adjacency[t].end()) {
        arc_adjacency[t][s] = 0;
      }
    } else {
      arc_adjacency[s][t] += cap;
    }
  };

  auto term_op = [&](bool is_source, int n, int cap) {
    g.nnode = std::max(g.nnode,n+1);
    g.terminal_capacities.resize(2 * (g.nnode), 0);
    if (is_source) {
      g.terminal_capacities[2 * n + 0] += cap;
    } else {
      g.terminal_capacities[2 * n + 1] += cap;
    }
  };

  _dimacs_implementation::read_dimacs_general(filename, arc_op, term_op);

  // process arcs after the fact
  g.narc = 0;
  for (const auto &[s, arc_list] : arc_adjacency) {
    for (const auto &[t, cap] : arc_list) {
      if (s < t) {
        g.arc_capacities.push_back(cap);
        g.arc_capacities.push_back(arc_adjacency[t][s]);
        g.arcs.push_back(s);
        g.arcs.push_back(t);
        g.narc++;
      }
    }
  }
  g.terminal_capacities.resize(2 * (g.nnode), 0); // ensure size
  g.arcs.shrink_to_fit();
  g.terminal_capacities.shrink_to_fit();
  g.arc_capacities.shrink_to_fit();
  return std::move(g);
}

} // namespace mcpd3
