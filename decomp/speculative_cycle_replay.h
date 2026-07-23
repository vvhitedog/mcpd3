// mcpd3 - minimum cut using a primal dual algorithm and dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <capacity.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mcpd3 {

struct BoundaryCycleSample {
  std::uint64_t state_hash = 0;
  std::vector<std::uint8_t> labels;
  std::vector<std::int8_t> diffs;
};

inline void record_lagrange_update_origin(const Lagrange &current_alpha,
                                          Lagrange &last_alpha,
                                          bool advance_last_alpha) {
  if (advance_last_alpha) {
    last_alpha = current_alpha;
  }
}

struct ExactBoundaryCycle {
  std::size_t period = 0;
  std::vector<BoundaryCycleSample> samples;

  const BoundaryCycleSample &nextSample(std::size_t offset) const {
    if (samples.empty() || period != samples.size()) {
      throw std::runtime_error("invalid exact boundary cycle");
    }
    return samples[offset % period];
  }
};

class ExactBoundaryCycleDetector {
public:
  ExactBoundaryCycleDetector(std::size_t max_period,
                             std::size_t minimum_repetitions)
      : max_period_(max_period),
        minimum_repetitions_(minimum_repetitions) {
    if (max_period_ < 2) {
      throw std::runtime_error(
          "speculative cycle maximum period must be at least two");
    }
    if (minimum_repetitions_ < 2) {
      throw std::runtime_error(
          "speculative cycle repetitions must be at least two");
    }
  }

  std::optional<ExactBoundaryCycle> observe(BoundaryCycleSample sample) {
    if (sample.labels.empty() || sample.diffs.empty()) {
      throw std::runtime_error(
          "speculative cycle samples must contain boundary state");
    }
    history_.push_back(std::move(sample));

    const std::size_t maximum_candidate =
        std::min(max_period_, history_.size() / minimum_repetitions_);
    for (std::size_t period = 2; period <= maximum_candidate; ++period) {
      if (!lastPeriodsMatch(period)) {
        continue;
      }
      ExactBoundaryCycle cycle;
      cycle.period = period;
      const std::size_t begin = history_.size() - period;
      cycle.samples.insert(cycle.samples.end(), history_.begin() + begin,
                           history_.end());
      trimHistory();
      return cycle;
    }

    trimHistory();
    return std::nullopt;
  }

  void clear() { history_.clear(); }

private:
  static bool sameState(const BoundaryCycleSample &lhs,
                        const BoundaryCycleSample &rhs) {
    return lhs.state_hash == rhs.state_hash && lhs.labels == rhs.labels;
  }

  bool lastPeriodsMatch(std::size_t period) const {
    const std::size_t period_begin = history_.size() - period;
    bool has_distinct_states = false;
    for (std::size_t offset = 1; offset < period; ++offset) {
      if (!sameState(history_[period_begin],
                     history_[period_begin + offset])) {
        has_distinct_states = true;
        break;
      }
    }
    if (!has_distinct_states) {
      return false;
    }

    for (std::size_t repetition = 1;
         repetition < minimum_repetitions_; ++repetition) {
      const std::size_t comparison_begin =
          period_begin - repetition * period;
      for (std::size_t offset = 0; offset < period; ++offset) {
        if (!sameState(history_[period_begin + offset],
                       history_[comparison_begin + offset])) {
          return false;
        }
      }
    }
    return true;
  }

  void trimHistory() {
    const std::size_t limit = max_period_ * minimum_repetitions_;
    if (history_.size() > limit) {
      history_.erase(history_.begin(),
                     history_.begin() + (history_.size() - limit));
    }
  }

  std::size_t max_period_;
  std::size_t minimum_repetitions_;
  std::vector<BoundaryCycleSample> history_;
};

inline bool accept_speculative_cycle_probe(
    const Objective &candidate_lower_bound,
    const Objective &baseline_lower_bound,
    const Objective &rollback_tolerance) {
  if (rollback_tolerance < 0) {
    throw std::runtime_error(
        "speculative cycle rollback tolerance must be non-negative");
  }
  if (candidate_lower_bound >= baseline_lower_bound) {
    return true;
  }
  return checked_subtract(
             baseline_lower_bound, candidate_lower_bound,
             "speculative cycle lower-bound difference overflow") <=
         rollback_tolerance;
}

} // namespace mcpd3
