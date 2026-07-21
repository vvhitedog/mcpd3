// mcpd3 - minimum cut using a primal dual algorithm and dual decomposition.
// Copyright (C) 2021 Matt Gara

#pragma once

namespace mcpd3 {

inline constexpr int kDefaultOptimizationScaleCount = 5;
inline constexpr int kDefaultMaxIterationCount = 10000;
inline constexpr long kDefaultInitialStepSize = 5000;
inline constexpr long kDefaultObjectiveScale = 500;
inline constexpr int kDefaultPatience = 10;
inline constexpr int kDefaultDisagreementPatience = 10;
inline constexpr bool kDefaultUseMomentum = true;
inline constexpr bool kDefaultEnableGroupStopping = false;
inline constexpr int kDefaultScaledEpsilonMaxStepSize = 12;
inline constexpr int kDefaultScaledEpsilonStrengthCap = 2;
inline constexpr int kDefaultMaxObjectiveScalePromotions = 4;
inline constexpr bool kDefaultRetryExhaustRegularizedScaleIterations = true;

} // namespace mcpd3
