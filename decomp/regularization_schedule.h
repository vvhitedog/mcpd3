// mcpd3 - minimum cut using a primal dual algorithm and dual decomposition.
// Copyright (C) 2021 Matt Gara

#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace mcpd3 {

enum class DualDecompositionRegularizationScheme {
  SCALED_EPSILON,
  DISAGREEMENT_PLATEAU_EPSILON,
  NONE
};

// Preserve the public worker name while making both execution paths consume
// exactly the same regularization policy type.
using PartitionWorkerRegularizationScheme =
    DualDecompositionRegularizationScheme;

template <typename Options>
long scaledEpsilonStrengthForStepSize(const Options &options,
                                      long step_size) {
  if (options.regularization_scheme !=
      DualDecompositionRegularizationScheme::SCALED_EPSILON) {
    return 0;
  }
  if (step_size > options.scaled_epsilon_max_step_size) {
    return 0;
  }
  return options.scaled_epsilon_strength_cap > 0
             ? std::min(step_size,
                        static_cast<long>(options.scaled_epsilon_strength_cap))
             : step_size;
}

class DisagreementPlateauRegularizationTracker {
public:
  explicit DisagreementPlateauRegularizationTracker(int patience)
      : patience_(patience) {
    if (patience_ <= 0) {
      throw std::runtime_error("disagreement patience must be positive");
    }
    reset();
  }

  void reset() {
    active_ = false;
    has_observation_ = false;
    has_regularization_pulse_ = false;
    best_disagreement_count_ = std::numeric_limits<long>::max();
    last_improvement_iteration_ = 0;
    last_regularization_pulse_iteration_ = 0;
  }

  bool observe(int iteration, long disagreement_count) {
    if (iteration < 0 || disagreement_count < 0) {
      throw std::runtime_error(
          "disagreement plateau observations must be non-negative");
    }
    if (!has_observation_ || disagreement_count < best_disagreement_count_) {
      has_observation_ = true;
      best_disagreement_count_ = disagreement_count;
      last_improvement_iteration_ = iteration;
      return false;
    }
    const int window_start =
        has_regularization_pulse_
            ? std::max(last_improvement_iteration_,
                       last_regularization_pulse_iteration_)
            : last_improvement_iteration_;
    if (iteration - window_start >= patience_) {
      active_ = true;
      has_regularization_pulse_ = true;
      last_regularization_pulse_iteration_ = iteration;
      return true;
    }
    return false;
  }

  bool active() const { return active_; }
  long bestDisagreementCount() const { return best_disagreement_count_; }
  int iterationsSinceImprovement(int iteration) const {
    return has_observation_ ? iteration - last_improvement_iteration_ : 0;
  }

private:
  int patience_;
  bool active_ = false;
  bool has_observation_ = false;
  bool has_regularization_pulse_ = false;
  long best_disagreement_count_ = std::numeric_limits<long>::max();
  int last_improvement_iteration_ = 0;
  int last_regularization_pulse_iteration_ = 0;
};

} // namespace mcpd3
