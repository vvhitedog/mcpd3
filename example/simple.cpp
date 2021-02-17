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
 * A simple example of computing a mincut using mcpd3 of the following graph,
 *
 *           SOURCE
 *          /       \
 *       10/         \5
 *        /      3    \
 *      node0 -----> node1
 *        |   <-----   |
 *        |      4     |
 *        \            /
 *        5\          /10
 *          \        /
 *             SINK
 *
 */

#include <iostream>
#include <primaldual/mcpd3.h>

int main(int argc, char *argv[]) {
  std::vector<int> arcs{0, 1};
  std::vector<int> arc_capacities{3, 4};
  std::vector<int> terminal_capacities{10, 5, 5, 10};
  mcpd3::PrimalDualMinCutSolver min_cut_solver(2, 1, std::move(arcs),
                                               std::move(arc_capacities),
                                               std::move(terminal_capacities));
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
