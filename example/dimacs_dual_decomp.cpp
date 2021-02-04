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
 * A example of a dual decomposition on a DIMACS format maxflow/mincut problem.
 */

#include <decomp/dualdecomp.h>
#include <graph/dimacs.h>
#include <iostream>

std::vector<int> basic_graph_partition(int npartition,
                                       const mcpd3::MinCutGraph &graph) {
  std::vector<int> partitions(graph.nnode);
  int nnode_in_each_partition = (graph.nnode + npartition - 1) / npartition;
  for (int i = 0; i < graph.nnode; ++i) {
    partitions[i] = i / nnode_in_each_partition;
  }
  return std::move(partitions);
}

int main(int argc, char *argv[]) {

  int npartition = 10;
  if (argc < 2) {
    std::cout << "usage: " << argv[0]
              << " DIMACS_MAXFLOW_FILE [NUM_PARTITIONS]\n";
    std::exit(EXIT_SUCCESS);
  }
  if (argc > 2) {
    npartition = std::atoi(argv[2]);
  }

  auto min_cut_graph_data = mcpd3::read_dimacs(argv[1]);
  //const int scale = 10000;
  //for ( auto &cap  :min_cut_graph_data.arc_capacities ) {
  //  cap *= scale;
  //}
  //for ( auto &cap  :min_cut_graph_data.terminal_capacities ) {
  //  cap *= scale;
  //}
  auto partitions = basic_graph_partition(npartition, min_cut_graph_data);
  mcpd3::DualDecomposition dual_decomp(npartition, std::move(partitions),
                                       std::move(min_cut_graph_data));
  //dual_decomp.runOptimizationStep(10000,10000,1);
  //dual_decomp.runOptimizationStep(10000,1000,2);
  //dual_decomp.runOptimizationStep(10000,100,3);
  //dual_decomp.runOptimizationStep(10000,10,4);
  //dual_decomp.runOptimizationStep(10000,1,5);
  //dual_decomp.runPrimalSolutionDecodingStep();
  //std::cout << "primal min cut value : " <<  dual_decomp.getPrimalMinCutValue() << "\n";


  dual_decomp.runOptimizationStep(10000,1,1);
  dual_decomp.scaleProblem<10>();
  dual_decomp.runOptimizationStep(10000,1,1);
  dual_decomp.scaleProblem<10>();
  dual_decomp.runOptimizationStep(10000,1,1);
  dual_decomp.scaleProblem<10>();
  dual_decomp.runOptimizationStep(10000,1,2);
  dual_decomp.scaleProblem<10>();
  dual_decomp.runOptimizationStep(10000,1,5);
  dual_decomp.scaleProblem<10>();
  dual_decomp.runOptimizationStep(10000,1,5);
  dual_decomp.runPrimalSolutionDecodingStep();
  std::cout << "primal min cut value : " <<  dual_decomp.getPrimalMinCutValue() << "\n";
  return EXIT_SUCCESS;
}
