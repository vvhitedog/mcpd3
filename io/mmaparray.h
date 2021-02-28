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

#include <fcntl.h>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <string>

namespace mcpd3 {

template <typename T> class MmapArray {
  static_assert(std::is_pod<T>::value, "MmapArray only supports POD types");

public:
  MmapArray(size_t num_elements, const std::string &filename)
      : num_elements_(num_elements), filename_(filename) {
    file_descriptor_ = open(filename_.c_str(), O_CREAT | O_RDWR, FILE_MODE);
    if (file_descriptor_ == -1) {
      throw std::runtime_error("failed to open file '" + filename_ + "'");
    }
    ftruncate(file_descriptor_, num_elements_ * sizeof(T));
    mem_map_ = reinterpret_cast<T *>(mmap(NULL, num_elements_ * sizeof(T),
                                          PROT_READ | PROT_WRITE, MAP_SHARED,
                                          file_descriptor_, 0));
  }

  ~MmapArray() {
    munmap(mem_map_, sizeof(T) * num_elements_);
    close(file_descriptor_);
  }

  T &operator[](size_t index) { return mem_map_[index]; }

  const T &operator[](size_t index) const { return mem_map_[index]; }

  T *data() { return mem_map_; }

  const T *data() const { return mem_map_; }

  T *begin() { return data(); }

  const T *begin() const { return data(); }

  T *end() { return data() + num_elements_; }

  const T *end() const { return data() + num_elements_; }

private:
  size_t num_elements_;
  std::string filename_;
  int file_descriptor_;

  T *mem_map_;

  const mode_t FILE_MODE = 0666;
};
} // namespace mcpd3
