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
#ifdef HAVE_METIS
#include <metis.h>
#endif

std::vector<int> basic_graph_partition(int npartition,
                                       const mcpd3::MinCutGraph &graph) {
  std::vector<int> partitions(graph.nnode);
  int nnode_in_each_partition = (graph.nnode + npartition - 1) / npartition;
  for (int i = 0; i < graph.nnode; ++i) {
    partitions[i] = i / nnode_in_each_partition;
  }
  return std::move(partitions);
}


#ifdef HAVE_METIS
std::vector<int> metis_partition(int npartition,
              const mcpd3::MinCutGraph &graph) {

  const auto &narc = graph.narc;
  const auto &arc = graph.arcs;
  const auto &nnode = graph.nnode;

  std::vector<std::vector<idx_t> > adj(nnode);

  // map an arc list to adjacency type that metis supports
  for (size_t aid = 0; aid < narc; ++aid) {
    int32_t s, t;
    s = arc[2 * aid + 0];
    t = arc[2 * aid + 1];
    adj[s].push_back(t);
    adj[t].push_back(s);
  }

  std::vector<idx_t> xadj(nnode + 1);
  std::vector<idx_t> adjv(narc * 2);

  size_t j = 0;
  size_t cumsum = 0;
  xadj[0] = cumsum;
  for (size_t inode = 0; inode < nnode; ++inode) {
    for (size_t i = 0; i < adj[inode].size(); ++i) {
      adjv[j++] = adj[inode][i];
    }
    cumsum += adj[inode].size();
    xadj[inode + 1] = cumsum;
  }

  std::vector<idx_t> part(nnode);
  idx_t ncon = 1;
  idx_t _nnode = nnode;
  idx_t _npart = npartition;
  idx_t objval;
  int ret = METIS_PartGraphKway(&_nnode, &ncon, &xadj[0], &adjv[0], nullptr,
                                nullptr, nullptr, &_npart, nullptr, nullptr,
                                nullptr, &objval, &part[0]);

  if (ret != METIS_OK) {
    fprintf(stderr, "Something failed while partitioning.\n");
    std::exit(EXIT_FAILURE);
  }

  std::vector<int> label(nnode,-1);
  std::copy(part.begin(), part.end(), label.begin());
  return std::move(label);
}
#endif

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

  //{ // run maxflow
  //  auto min_cut_graph_data_copy = min_cut_graph_data;
  //  mcpd3::PrimalDualMinCutSolver solver(std::move(min_cut_graph_data_copy));
  //  auto microseconds = mcpd3::time_lambda([&]{
  //      solver.maxflow(); // compute maxflow
  //      });
  //  std::cout << " maxflow : " << microseconds.count() << "ms\n";
  //}
#ifdef HAVE_METIS
  auto partitions = metis_partition(npartition, min_cut_graph_data);
#else
  auto partitions = basic_graph_partition(npartition, min_cut_graph_data);
#endif

  mcpd3::DualDecomposition dual_decomp(npartition, std::move(partitions),
                                       std::move(min_cut_graph_data));
  auto microseconds = mcpd3::time_lambda( [&] {
  const int scaling_factor = 10;
  const int num_optimization_scales = 6;
  for ( int iscale = 0; iscale < num_optimization_scales; ++iscale ) {
    auto status = dual_decomp.runOptimizationScale(10000,1,1,true);
    //dual_decomp.runPrimalSolutionDecodingStep(true);
    if ( status == mcpd3::DualDecomposition::OPTIMAL ) {
      break;
    }
    dual_decomp.scaleProblem<scaling_factor>();
  }
  });
    std::cout << " full loop time : " << microseconds.count() << "ms\n";

  //dual_decomp.runOptimizationScale(10000,1,1,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,1,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,2,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,5,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,5,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,5,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  //dual_decomp.scaleProblem<scaling_factor>();
  //dual_decomp.runOptimizationScale(10000,1,5,true);
  //dual_decomp.runPrimalSolutionDecodingStep(true);
  dual_decomp.runPrimalSolutionDecodingStep();
  std::cout << "primal min cut value : " <<  dual_decomp.getPrimalMinCutValue() << "\n";
  std::cout << " total solve loop time: " << dual_decomp.getTotalSolveLoopTime() << "\n";
  return EXIT_SUCCESS;
}
