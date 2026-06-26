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

#include <graph/mcgraph.h>
#include <vector>
#ifdef HAVE_METIS
#include <metis.h>
#endif

namespace mcpd3 {

std::vector<int> basic_graph_partition(int npartition, int narc, int nnode,
                                       const std::vector<int> &arc) {
  std::vector<int> partitions(nnode);
  int nnode_in_each_partition = (nnode + npartition - 1) / npartition;
  for (int i = 0; i < nnode; ++i) {
    partitions[i] = i / nnode_in_each_partition;
  }
  return std::move(partitions);
}

std::vector<int> basic_graph_partition(int npartition,
                                       const mcpd3::MinCutGraph &graph) {
  return basic_graph_partition(npartition, graph.narc, graph.nnode, graph.arcs);
}

#ifdef HAVE_METIS

std::vector<int> metis_partition(int npartition, int narc, int nnode,
                                 const std::vector<int> &arc) {

  std::vector<idx_t> xadj(nnode + 1);
  for (size_t aid = 0; aid < static_cast<size_t>(narc); ++aid) {
    const int32_t s = arc[2 * aid + 0];
    const int32_t t = arc[2 * aid + 1];
    ++xadj[s + 1];
    ++xadj[t + 1];
  }
  for (int inode = 0; inode < nnode; ++inode) {
    xadj[inode + 1] += xadj[inode];
  }

  std::vector<idx_t> adjv(static_cast<size_t>(xadj[nnode]));
  std::vector<idx_t> cursor = xadj;
  for (size_t aid = 0; aid < static_cast<size_t>(narc); ++aid) {
    const int32_t s = arc[2 * aid + 0];
    const int32_t t = arc[2 * aid + 1];
    adjv[cursor[s]++] = t;
    adjv[cursor[t]++] = s;
  }
  cursor.clear();
  cursor.shrink_to_fit();

  std::vector<idx_t> part(nnode);
  idx_t ncon = 1;
  idx_t _nnode = nnode;
  idx_t _npart = npartition;
  idx_t objval;
  int ret = METIS_PartGraphKway(&_nnode, &ncon, &xadj[0], &adjv[0], nullptr,
                                nullptr, nullptr, &_npart, nullptr, nullptr,
                                nullptr, &objval, &part[0]);

  if (ret != METIS_OK) {
    throw std::runtime_error("Something failed while partitioning.\n");
  }

  std::vector<int> label(nnode, -1);
  std::copy(part.begin(), part.end(), label.begin());
  return std::move(label);
}

std::vector<int> metis_partition(int npartition,
                                 const mcpd3::MinCutGraph &graph) {
  const auto &narc = graph.narc;
  const auto &arc = graph.arcs;
  const auto &nnode = graph.nnode;
  return metis_partition(npartition, narc, nnode, arc);
}
#endif
} // namespace mcpd3
