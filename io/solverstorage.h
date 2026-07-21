// mcpd3 - minimum cut using a primal dual algorithm and dual decomposition.
// Copyright (C) 2021 Matt Gara

#pragma once

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <vector>

namespace mcpd3 {

enum class SolverStorageMode {
  RESIDENT,
  FILE_BACKED_MMAP,
  ANONYMOUS_MMAP,
};

struct SolverStorageOptions {
  SolverStorageMode mode = SolverStorageMode::RESIDENT;
  std::string directory;
  std::string mmap_advice;
};

template <typename T>
inline constexpr bool solver_storage_mmap_compatible_v =
    std::is_trivially_copyable_v<T>;

inline void applySolverMmapAdvice(void *address, std::size_t bytes,
                                  const SolverStorageOptions &options,
                                  const std::string &kind) {
  if (address == nullptr || bytes == 0 || options.mmap_advice.empty() ||
      options.mmap_advice == "none") {
    return;
  }
  int result = 0;
  if (options.mmap_advice == "normal") {
    result = madvise(address, bytes, MADV_NORMAL);
  } else if (options.mmap_advice == "random") {
    result = madvise(address, bytes, MADV_RANDOM);
  } else if (options.mmap_advice == "sequential") {
    result = madvise(address, bytes, MADV_SEQUENTIAL);
  } else if (options.mmap_advice == "willneed" ||
             options.mmap_advice == "populate") {
    result = madvise(address, bytes, MADV_WILLNEED);
  } else if (options.mmap_advice == "lock") {
    result = mlock(address, bytes);
  } else if (options.mmap_advice == "lock_onfault") {
#ifdef MLOCK_ONFAULT
    result = mlock2(address, bytes, MLOCK_ONFAULT);
#else
    throw std::invalid_argument(
        "solver mmap advice lock_onfault is unavailable on this platform");
#endif
  } else if (options.mmap_advice == "dontdump") {
#ifdef MADV_DONTDUMP
    result = madvise(address, bytes, MADV_DONTDUMP);
#endif
  } else {
    throw std::invalid_argument("unknown solver mmap advice '" +
                                options.mmap_advice + "'");
  }
  if (result != 0) {
    throw std::runtime_error("madvise failed for solver " + kind + ": " +
                             std::strerror(errno));
  }
}

template <typename T> class SolverArray {
public:
  SolverArray() = default;

  SolverArray(std::size_t count, const T &value,
              const SolverStorageOptions &options, const std::string &kind) {
    initialize(count, options, kind);
    std::fill(begin(), end(), value);
  }

  SolverArray(std::size_t count, const SolverStorageOptions &options,
              const std::string &kind) {
    initialize(count, options, kind);
  }

  SolverArray(std::vector<T> values, const SolverStorageOptions &options,
              const std::string &kind) {
    if (options.mode == SolverStorageMode::RESIDENT) {
      resident_ = std::move(values);
      size_ = resident_.size();
      return;
    }
    initialize(values.size(), options, kind);
    std::copy(values.begin(), values.end(), begin());
  }

  SolverArray(const SolverArray &) = delete;
  SolverArray &operator=(const SolverArray &) = delete;

  SolverArray(SolverArray &&other) noexcept { moveFrom(std::move(other)); }

  SolverArray &operator=(SolverArray &&other) noexcept {
    if (this != &other) {
      release();
      moveFrom(std::move(other));
    }
    return *this;
  }

  ~SolverArray() { release(); }

  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  T *data() { return mapped_ != nullptr ? mapped_ : resident_.data(); }
  const T *data() const {
    return mapped_ != nullptr ? mapped_ : resident_.data();
  }
  T *begin() { return data(); }
  const T *begin() const { return data(); }
  T *end() { return size_ == 0 ? data() : data() + size_; }
  const T *end() const { return size_ == 0 ? data() : data() + size_; }
  T &operator[](std::size_t index) { return data()[index]; }
  const T &operator[](std::size_t index) const { return data()[index]; }

  bool isFileBacked() const {
    return mapped_ != nullptr &&
           mode_ == SolverStorageMode::FILE_BACKED_MMAP;
  }
  std::size_t fileBackedBytes() const {
    return isFileBacked() ? mapped_bytes_ : 0;
  }

  std::vector<T> copyToVector() const {
    return std::vector<T>(begin(), end());
  }

  void replace(const std::vector<T> &values) {
    requireSameSize(values.size());
    std::copy(values.begin(), values.end(), begin());
  }

  void replace(std::vector<T> &&values) {
    requireSameSize(values.size());
    if (mode_ == SolverStorageMode::RESIDENT) {
      resident_ = std::move(values);
      size_ = resident_.size();
      return;
    }
    std::move(values.begin(), values.end(), begin());
  }

  template <typename Container> void replaceFrom(const Container &values) {
    requireSameSize(values.size());
    std::copy(values.begin(), values.end(), begin());
  }

  bool equals(const std::vector<T> &values) const {
    return values.size() == size_ &&
           std::equal(begin(), end(), values.begin());
  }

  bool equals(const SolverArray &values) const {
    return values.size() == size_ &&
           std::equal(begin(), end(), values.begin());
  }

  SolverArray clone(const SolverStorageOptions &options,
                    const std::string &kind) const {
    SolverArray copy(size_, options, kind);
    std::copy(begin(), end(), copy.begin());
    return copy;
  }

private:
  void initialize(std::size_t count, const SolverStorageOptions &options,
                  const std::string &kind) {
    size_ = count;
    mode_ = options.mode;
    if (mode_ == SolverStorageMode::RESIDENT) {
      resident_.resize(count);
      return;
    }
    if constexpr (!solver_storage_mmap_compatible_v<T>) {
      throw std::invalid_argument(
          "file-backed solver arrays require trivially copyable values");
    } else {
      if (count == 0) {
        return;
      }
      mapped_bytes_ = count * sizeof(T);
      int flags = MAP_PRIVATE | MAP_ANONYMOUS;
      if (mode_ == SolverStorageMode::FILE_BACKED_MMAP) {
        if (options.directory.empty()) {
          throw std::invalid_argument(
              "file-backed solver storage requires a directory");
        }
        std::string pattern = options.directory + "/mcpd3_" + kind + "_XXXXXX";
        fd_ = mkstemp(pattern.data());
        if (fd_ == -1) {
          throw std::runtime_error("failed to create solver mmap file '" +
                                   pattern + "': " + std::strerror(errno));
        }
        unlink(pattern.c_str());
        if (ftruncate(fd_, static_cast<off_t>(mapped_bytes_)) != 0) {
          const std::string message =
              "failed to size solver mmap file: " +
              std::string(std::strerror(errno));
          close(fd_);
          fd_ = -1;
          throw std::runtime_error(message);
        }
        flags = MAP_SHARED;
      }
      void *address = mmap(nullptr, mapped_bytes_, PROT_READ | PROT_WRITE,
                           flags, fd_, 0);
      if (address == MAP_FAILED) {
        const std::string message =
            "failed to mmap solver " + kind + ": " + std::strerror(errno);
        if (fd_ != -1) {
          close(fd_);
          fd_ = -1;
        }
        throw std::runtime_error(message);
      }
      mapped_ = static_cast<T *>(address);
      try {
        applySolverMmapAdvice(mapped_, mapped_bytes_, options, kind);
      } catch (...) {
        munmap(mapped_, mapped_bytes_);
        mapped_ = nullptr;
        if (fd_ != -1) {
          close(fd_);
          fd_ = -1;
        }
        throw;
      }
    }
  }

  void requireSameSize(std::size_t size) const {
    if (size != size_) {
      throw std::invalid_argument("solver array replacement size mismatch");
    }
  }

  void release() noexcept {
    if (mapped_ != nullptr) {
      munmap(mapped_, mapped_bytes_);
      mapped_ = nullptr;
    }
    if (fd_ != -1) {
      close(fd_);
      fd_ = -1;
    }
    resident_.clear();
    size_ = 0;
    mapped_bytes_ = 0;
  }

  void moveFrom(SolverArray &&other) noexcept {
    resident_ = std::move(other.resident_);
    mapped_ = other.mapped_;
    size_ = other.size_;
    mapped_bytes_ = other.mapped_bytes_;
    fd_ = other.fd_;
    mode_ = other.mode_;
    other.mapped_ = nullptr;
    other.size_ = 0;
    other.mapped_bytes_ = 0;
    other.fd_ = -1;
  }

  std::vector<T> resident_;
  T *mapped_ = nullptr;
  std::size_t size_ = 0;
  std::size_t mapped_bytes_ = 0;
  int fd_ = -1;
  SolverStorageMode mode_ = SolverStorageMode::RESIDENT;
};

} // namespace mcpd3
