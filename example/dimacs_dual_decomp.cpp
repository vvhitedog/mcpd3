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

#include <io/memory.h>

#include <decomp/dualdecomp.h>
#include <graph/dimacs.h>
#include <iostream>

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
  mcpd3::MinCutGraph min_cut_graph_data;
  auto read_graph_mem = mcpd3::get_resident_memory_usage(
      [&] { min_cut_graph_data = mcpd3::read_dimacs(argv[1]); });
  std::cout << "read graph mem usage: " << read_graph_mem.usage_in_gb << "GB\n";

  mcpd3::DualDecomposition dual_decomp(npartition, 
                                       std::move(min_cut_graph_data));
  auto microseconds = mcpd3::time_lambda([&] {
    const int scaling_factor = 10;
    const int num_optimization_scales = 6;
    for (int iscale = 0; iscale < num_optimization_scales; ++iscale) {
      auto status = dual_decomp.runOptimizationScale(10000, 1, 1, true);
      // dual_decomp.runPrimalSolutionDecodingStep(true);
      if (status == mcpd3::DualDecomposition::OPTIMAL) {
        break;
      }
      dual_decomp.scaleProblem<scaling_factor>();
    }
  });
  std::cout << " full loop time : " << microseconds.count() << "ms\n";

  // dual_decomp.runOptimizationScale(10000,1,1,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,1,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,2,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,5,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,5,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,5,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.scaleProblem<scaling_factor>();
  // dual_decomp.runOptimizationScale(10000,1,5,true);
  // dual_decomp.runPrimalSolutionDecodingStep(true);
  // dual_decomp.runPrimalSolutionDecodingStep();
  // std::cout << "primal min cut value : " <<
  // dual_decomp.getPrimalMinCutValue() << "\n";
  std::cout << " total solve loop time: " << dual_decomp.getTotalSolveLoopTime()
            << "\n";
  return EXIT_SUCCESS;
}
