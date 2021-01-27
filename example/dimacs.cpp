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

/**
 * A example of computing a mincut using mcpd3 on dimacs format graphs.
 */

#include <iostream>
#include <primaldual/mcpd3.h>
#include <unordered_map>

struct RawGraphData {
  int nnode, narc;
  std::vector<int> arcs;
  std::vector<int> arc_capacities;
  std::vector<int> terminal_capacities;
};

RawGraphData read_dimacs(const std::string &filename) {
  RawGraphData g;
  std::unordered_map<int, int> remapped_indices;
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
      g.terminal_capacities.resize(2 * (n - 2), 0);
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
      if (s != source && remapped_indices.find(s) == remapped_indices.end()) {
        remapped_indices.insert({s,remapped_indices.size()});
      }
      if (t != sink && remapped_indices.find(t) == remapped_indices.end()) {
        remapped_indices.insert({t,remapped_indices.size()});
      }
      if (s != source && t != sink) { // handle non terminal arc
        _s = remapped_indices[s];
        _t = remapped_indices[t];
        if (arc_adjacency[_s].find(_t) == arc_adjacency[_s].end()) {
          arc_adjacency[_s][_t] = cap;
          if (arc_adjacency[_t].find(_s) == arc_adjacency[_t].end()) {
            arc_adjacency[_t][_s] = 0;
          }
        } else {
          arc_adjacency[_s][_t] += cap;
        }
      } else { // handle terminal arc
        if (s == source) {
          _t = remapped_indices[t];
          g.terminal_capacities[2 * _t + 0] += cap;
        } else if (t == sink) {
          _s = remapped_indices[s];
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
  assert(remapped_indices.find(source) == remapped_indices.end());
  assert(remapped_indices.find(sink) == remapped_indices.end());
  g.nnode = remapped_indices.size();
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

int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cout << "usage: " << argv[0] << " DIMACS_MAXFLOW_FILE\n";
    std::exit(EXIT_SUCCESS);
  }

  auto raw_graph_data = read_dimacs(argv[1]);

  mcpd3::PrimalDualMinCutSolver min_cut_solver(
      raw_graph_data.nnode, raw_graph_data.narc, std::move(raw_graph_data.arcs),
      std::move(raw_graph_data.arc_capacities),
      std::move(raw_graph_data.terminal_capacities));
  std::cout << " ========= direct maxflow computation =========\n";
  auto maxflow_value = min_cut_solver.maxflow();
  std::cout << " maxflow value: " << maxflow_value << "\n";
  std::cout << " ========= primal dual maxflow computation ====\n";
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  std::cout << " ========= primal dual maxflow re-computation==\n";
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  std::cout << " ==============================================\n";
  std::cout << " ========= all values above should match ======\n";
  std::cout << " ==============================================\n";
  return EXIT_SUCCESS;
}
