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

#include <unordered_map>
#include <string>
#include <cassert>

#include <graph/mcgraph.h>

namespace mcpd3 {

int remap_index(int original_index, int source_index, int sink_index ) {
  if ( original_index > source_index ) {
    --original_index;
  }
  if ( original_index > sink_index ) {
    --original_index;
  }
  return --original_index; // indices in DIMACS start at 1 but we wish for them to start at 0
}

MinCutGraph read_dimacs(const std::string &filename) {
  MinCutGraph g;
  std::unordered_map<int, std::unordered_map<int, int>> arc_adjacency;
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
      g.nnode = n-2;
      g.terminal_capacities.resize(2*(g.nnode),0);
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
        _s = remap_index(s,source,sink);
        _t = remap_index(t,source,sink);
        if (arc_adjacency[_s].find(_t) == arc_adjacency[_s].end()) {
          arc_adjacency[_s][_t] = cap;
          if (arc_adjacency[_t].find(_s) == arc_adjacency[_t].end()) {
            arc_adjacency[_t][_s] = 0;
          }
        } else {
          arc_adjacency[_s][_t] = cap;
        }
      } else { // handle terminal arc
        if (s == source) {
          _t = remap_index(t,source,sink);
          g.terminal_capacities[2 * _t + 0] += cap;
        } else if (t == sink) {
          _s = remap_index(s,source,sink);
          g.terminal_capacities[2 * _s + 1] += cap;
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
  fclose(stream);
  return std::move(g);

}

}
