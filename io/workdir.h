// mcpd3 - minimum cut using a primal dual algorithm and the dual decomposition.
// Copyright (C) 2021 Matt Gara
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.

#pragma once

#include <filesystem>
#include <string>

namespace mcpd3 {

inline void create_work_directory(const std::string &path) {
  std::filesystem::remove_all(path);
  std::filesystem::create_directories(path);
}

} // namespace mcpd3
