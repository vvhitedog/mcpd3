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
#include <graph/dimacs.h>


int main(int argc, char *argv[]) {

  if (argc < 2) {
    std::cout << "usage: " << argv[0] << " DIMACS_MAXFLOW_FILE\n";
    std::exit(EXIT_SUCCESS);
  }

  auto min_cut_graph_data = mcpd3::read_dimacs(argv[1]);
  mcpd3::PrimalDualMinCutSolver min_cut_solver(std::move(min_cut_graph_data));
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
