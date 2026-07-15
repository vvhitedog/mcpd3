// mcpd3 - minimum cut using a primal dual algorithm and dual decomposition.
// Copyright (C) 2021 Matt Gara

#pragma once

#include <algorithm>
#include <deque>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace mcpd3 {

inline constexpr int kInfiniteHaloDepth = -1;

struct HaloPartitionLayout {
  long objective_multiplier = 1;
  std::vector<std::vector<int>> node_partitions;
  std::vector<std::vector<int>> arc_partitions;
};

inline long checkedHaloLcm(long lhs, int rhs) {
  if (lhs <= 0 || rhs <= 0) {
    throw std::runtime_error("halo multiplicities must be positive");
  }
  const long divisor = std::gcd(lhs, static_cast<long>(rhs));
  const long quotient = lhs / divisor;
  if (quotient > std::numeric_limits<long>::max() / rhs) {
    throw std::overflow_error("halo objective multiplier overflow");
  }
  return quotient * rhs;
}

inline void insertHaloMembership(std::vector<int> *memberships,
                                 int partition) {
  const auto position =
      std::lower_bound(memberships->begin(), memberships->end(), partition);
  if (position == memberships->end() || *position != partition) {
    memberships->insert(position, partition);
  }
}

inline HaloPartitionLayout buildHaloPartitionLayout(
    int partition_count, int node_count, const std::vector<int> &arcs,
    const std::vector<int> &partition_labels, int halo_depth) {
  if (partition_count <= 0) {
    throw std::runtime_error("halo partition count must be positive");
  }
  if (node_count < 0) {
    throw std::runtime_error("halo node count must be non-negative");
  }
  if (halo_depth != kInfiniteHaloDepth && halo_depth < 1) {
    throw std::runtime_error(
        "halo depth must be positive or kInfiniteHaloDepth");
  }
  if (partition_labels.size() != static_cast<size_t>(node_count)) {
    throw std::runtime_error(
        "halo partition label count must match node count");
  }
  if (arcs.size() % 2 != 0) {
    throw std::runtime_error("halo arcs must contain endpoint pairs");
  }
  for (const int partition : partition_labels) {
    if (partition < 0 || partition >= partition_count) {
      throw std::runtime_error("halo partition label is out of range");
    }
  }
  for (const int node : arcs) {
    if (node < 0 || node >= node_count) {
      throw std::runtime_error("halo arc endpoint is out of range");
    }
  }

  const size_t arc_count = arcs.size() / 2;
  HaloPartitionLayout layout;
  layout.node_partitions.resize(static_cast<size_t>(node_count));
  layout.arc_partitions.resize(arc_count);

  // Preserve the compact historical mcpd3-n representation exactly at h1.
  // Each canonically oriented edge has one owner, and only endpoints needed
  // by those owned edges become local node copies.
  if (halo_depth == 1) {
    for (size_t arc = 0; arc < arc_count; ++arc) {
      int source = arcs[2 * arc];
      int target = arcs[2 * arc + 1];
      if (source > target) {
        std::swap(source, target);
      }
      const int owner = partition_labels[static_cast<size_t>(source)];
      layout.arc_partitions[arc].push_back(owner);
      insertHaloMembership(
          &layout.node_partitions[static_cast<size_t>(source)], owner);
      insertHaloMembership(
          &layout.node_partitions[static_cast<size_t>(target)], owner);
      insertHaloMembership(
          &layout.node_partitions[static_cast<size_t>(target)],
          partition_labels[static_cast<size_t>(target)]);
    }
    return layout;
  }

  if (halo_depth == kInfiniteHaloDepth) {
    for (auto &memberships : layout.node_partitions) {
      memberships.resize(static_cast<size_t>(partition_count));
      std::iota(memberships.begin(), memberships.end(), 0);
    }
    for (auto &memberships : layout.arc_partitions) {
      memberships.resize(static_cast<size_t>(partition_count));
      std::iota(memberships.begin(), memberships.end(), 0);
    }
    layout.objective_multiplier = partition_count;
    return layout;
  }

  std::vector<std::vector<int>> core_nodes(
      static_cast<size_t>(partition_count));
  for (int node = 0; node < node_count; ++node) {
    core_nodes[static_cast<size_t>(partition_labels[node])].push_back(node);
  }

  std::vector<int> offsets(static_cast<size_t>(node_count) + 1, 0);
  for (size_t arc = 0; arc < arc_count; ++arc) {
    ++offsets[static_cast<size_t>(arcs[2 * arc]) + 1];
    ++offsets[static_cast<size_t>(arcs[2 * arc + 1]) + 1];
  }
  for (int node = 0; node < node_count; ++node) {
    offsets[static_cast<size_t>(node) + 1] +=
        offsets[static_cast<size_t>(node)];
  }
  std::vector<int> adjacency(static_cast<size_t>(offsets.back()));
  std::vector<int> cursor = offsets;
  for (size_t arc = 0; arc < arc_count; ++arc) {
    const int source = arcs[2 * arc];
    const int target = arcs[2 * arc + 1];
    adjacency[static_cast<size_t>(cursor[source]++)] = target;
    adjacency[static_cast<size_t>(cursor[target]++)] = source;
  }

  std::vector<int> distance(static_cast<size_t>(node_count), -1);
  std::vector<int> visited;
  std::deque<int> queue;
  for (int partition = 0; partition < partition_count; ++partition) {
    queue.clear();
    visited.clear();
    for (const int node : core_nodes[static_cast<size_t>(partition)]) {
      distance[static_cast<size_t>(node)] = 0;
      visited.push_back(node);
      queue.push_back(node);
    }
    while (!queue.empty()) {
      const int node = queue.front();
      queue.pop_front();
      insertHaloMembership(
          &layout.node_partitions[static_cast<size_t>(node)], partition);
      if (distance[static_cast<size_t>(node)] == halo_depth) {
        continue;
      }
      for (int index = offsets[static_cast<size_t>(node)];
           index < offsets[static_cast<size_t>(node) + 1]; ++index) {
        const int neighbor = adjacency[static_cast<size_t>(index)];
        if (distance[static_cast<size_t>(neighbor)] >= 0) {
          continue;
        }
        distance[static_cast<size_t>(neighbor)] =
            distance[static_cast<size_t>(node)] + 1;
        visited.push_back(neighbor);
        queue.push_back(neighbor);
      }
    }
    for (const int node : visited) {
      distance[static_cast<size_t>(node)] = -1;
    }
  }

  for (size_t arc = 0; arc < arc_count; ++arc) {
    const auto &source_memberships =
        layout.node_partitions[static_cast<size_t>(arcs[2 * arc])];
    const auto &target_memberships =
        layout.node_partitions[static_cast<size_t>(arcs[2 * arc + 1])];
    std::set_intersection(
        source_memberships.begin(), source_memberships.end(),
        target_memberships.begin(), target_memberships.end(),
        std::back_inserter(layout.arc_partitions[arc]));
    if (layout.arc_partitions[arc].empty()) {
      throw std::runtime_error("halo construction omitted an original edge");
    }
  }

  for (const auto &memberships : layout.node_partitions) {
    if (memberships.empty()) {
      throw std::runtime_error("halo construction omitted an original node");
    }
    layout.objective_multiplier = checkedHaloLcm(
        layout.objective_multiplier, static_cast<int>(memberships.size()));
  }
  for (const auto &memberships : layout.arc_partitions) {
    layout.objective_multiplier = checkedHaloLcm(
        layout.objective_multiplier, static_cast<int>(memberships.size()));
  }
  return layout;
}

} // namespace mcpd3
