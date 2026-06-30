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

#include <chrono>

namespace mcpd3 {

/**
 * @brief Time a lambda function.
 *
 * @param lambda - the function to execute and time
 *
 * @return the number of microseconds elapsed while executing lambda
 */
template <typename Lambda>
std::chrono::microseconds time_lambda(Lambda lambda) {
  auto start_time = std::chrono::high_resolution_clock::now();
  lambda();
  auto end_time = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end_time -
                                                               start_time);
}

}
