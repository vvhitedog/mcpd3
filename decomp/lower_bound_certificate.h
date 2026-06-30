#pragma once

#include <limits>
#include <stdexcept>

namespace mcpd3 {

inline long checkedAddObjectiveRaw(long lhs, long rhs, const char *context) {
  if (rhs > 0 && lhs > std::numeric_limits<long>::max() - rhs) {
    throw std::overflow_error(context);
  }
  if (rhs < 0 && lhs < std::numeric_limits<long>::min() - rhs) {
    throw std::overflow_error(context);
  }
  return lhs + rhs;
}

inline long checkedSubtractObjectiveRaw(long lhs, long rhs,
                                        const char *context) {
  if (rhs > 0 && lhs < std::numeric_limits<long>::min() + rhs) {
    throw std::overflow_error(context);
  }
  if (rhs < 0 && lhs > std::numeric_limits<long>::max() + rhs) {
    throw std::overflow_error(context);
  }
  return lhs - rhs;
}

inline long regularizedObjectiveRaw(long original_objective_raw,
                                    long regularization_contribution_raw) {
  return checkedAddObjectiveRaw(original_objective_raw,
                                regularization_contribution_raw,
                                "regularized objective overflow");
}

inline long certifiedOriginalLowerBoundRaw(
    long original_objective_raw, long regularization_contribution_raw,
    long regularization_budget_raw) {
  return checkedSubtractObjectiveRaw(
      regularizedObjectiveRaw(original_objective_raw,
                              regularization_contribution_raw),
      regularization_budget_raw, "certified lower bound overflow");
}

} // namespace mcpd3
