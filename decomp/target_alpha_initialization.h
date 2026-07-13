#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace mcpd3 {

inline std::vector<double> targetAlphaUnaryProposal(
    long target_value, long current_value,
    const std::vector<int> &target_labels,
    const std::vector<int> &current_labels,
    const std::vector<unsigned char> &eligible, double damping) {
  if (target_labels.size() != current_labels.size() ||
      target_labels.size() != eligible.size()) {
    throw std::invalid_argument(
        "target alpha proposal vectors must have equal sizes");
  }
  if (!std::isfinite(damping) || damping < 0.0) {
    throw std::invalid_argument(
        "target alpha proposal damping must be finite and nonnegative");
  }

  std::vector<double> proposal(target_labels.size(), 0.0);
  if (damping == 0.0 || target_value <= current_value) {
    return proposal;
  }
  const long double objective_deficit =
      static_cast<long double>(target_value) -
      static_cast<long double>(current_value);

  long difference_count = 0;
  for (size_t index = 0; index < target_labels.size(); ++index) {
    if ((target_labels[index] != 0 && target_labels[index] != 1) ||
        (current_labels[index] != 0 && current_labels[index] != 1)) {
      throw std::invalid_argument(
          "target alpha proposal labels must be binary");
    }
    if (eligible[index] && target_labels[index] != current_labels[index]) {
      ++difference_count;
    }
  }
  if (difference_count == 0) {
    return proposal;
  }

  const double per_difference = static_cast<double>(
      static_cast<long double>(damping) * objective_deficit /
      static_cast<long double>(difference_count));
  if (!std::isfinite(per_difference)) {
    throw std::overflow_error("target alpha proposal overflow");
  }
  for (size_t index = 0; index < target_labels.size(); ++index) {
    if (!eligible[index]) {
      continue;
    }
    const int difference = target_labels[index] - current_labels[index];
    proposal[index] = -per_difference * static_cast<double>(difference);
  }
  return proposal;
}

} // namespace mcpd3
