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

std::vector<int> basic_graph_partition(int npartition,
                                 int narc,
                                 int nnode,
                                 const std::vector<int> &arc ) {
  std::vector<int> partitions(nnode);
  int nnode_in_each_partition = (nnode + npartition - 1) / npartition;
  for (int i = 0; i < nnode; ++i) {
    partitions[i] = i / nnode_in_each_partition;
  }
  return std::move(partitions);
}

std::vector<int> basic_graph_partition(int npartition,
                                       const mcpd3::MinCutGraph &graph) {
  return basic_graph_partition(npartition,graph.narc,graph.nnode,graph.arcs);
}

#ifdef HAVE_METIS

std::vector<int> metis_partition(int npartition,
                                 int narc,
                                 int nnode,
                                 const std::vector<int> &arc ) {


  std::vector<std::vector<idx_t>> adj(nnode);

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
    throw std::runtime_error( "Something failed while partitioning.\n");
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
  return metis_partition(npartition,narc,nnode,arc);
}
#endif
}
