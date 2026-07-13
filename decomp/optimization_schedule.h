// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

namespace mcpd3 {

inline long nextOptimizationScheduleValue(long current) {
  if (current <= 1) {
    return 0;
  }
  const long next = current / 10;
  return next >= 1 ? next : 1;
}

} // namespace mcpd3
