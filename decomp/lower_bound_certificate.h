#pragma once

#include <capacity.h>

namespace mcpd3 {

inline Objective checkedAddObjectiveRaw(const Objective &lhs,
                                        const Objective &rhs,
                                        const char *context) {
  return checked_add(lhs, rhs, context);
}

inline Objective checkedSubtractObjectiveRaw(const Objective &lhs,
                                             const Objective &rhs,
                                             const char *context) {
  return checked_subtract(lhs, rhs, context);
}

inline Objective regularizedObjectiveRaw(
    const Objective &original_objective_raw,
    const Objective &regularization_contribution_raw) {
  return checkedAddObjectiveRaw(original_objective_raw,
                                regularization_contribution_raw,
                                "regularized objective overflow");
}

inline Objective certifiedOriginalLowerBoundRaw(
    const Objective &original_objective_raw,
    const Objective &regularization_contribution_raw,
    const Objective &regularization_budget_raw) {
  return checkedSubtractObjectiveRaw(
      regularizedObjectiveRaw(original_objective_raw,
                              regularization_contribution_raw),
      regularization_budget_raw, "certified lower bound overflow");
}

} // namespace mcpd3
