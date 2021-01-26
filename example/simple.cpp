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

#include <primaldual/mcpd3.h>
#include <iostream>

int main(int argc, char *argv[]) {
  {
  std::vector<int> arcs{0, 1};
  std::vector<int> arc_capacities{3, 4};
  std::vector<int> terminal_capacities{10, 5, 5, 10};

  mcpd3::PrimalDualMinCutSolver min_cut_solver(2, 1, std::move(arcs),
                                std::move(arc_capacities),
                                std::move(terminal_capacities));
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  }

  {
  std::vector<int> arcs{0, 1,
  0, 2,
  2, 3,
  1, 3
  };
  std::vector<int> arc_capacities{500, 0,
  200, 0,
  100, 0,
  0, 800
  };
  std::vector<int> terminal_capacities{100, 0,
    0, 70,
    80, 0,
    0, 100
  };

  mcpd3::PrimalDualMinCutSolver min_cut_solver(4, 4, std::move(arcs),
                                std::move(arc_capacities),
                                std::move(terminal_capacities));
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  min_cut_solver.terminal_capacities()[0] += 100;
  min_cut_solver.terminal_capacities()[3] += 100;
  min_cut_solver.terminal_capacities()[4] += 100;
  min_cut_solver.terminal_capacities()[7] += 100;
  min_cut_solver.solve();
  std::cout << " min cut value: " << min_cut_solver.getMinCutValue() << "\n";
  std::cout << " maxflow value: " << min_cut_solver.getMaxFlowValue() << "\n";
  }
  return EXIT_SUCCESS;
}
