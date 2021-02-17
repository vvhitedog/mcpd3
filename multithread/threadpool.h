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

#include <condition_variable>
#include <future>
#include <list>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace mcpd3 {

template <typename F> class ThreadPool {

public:
  ThreadPool(const size_t thread_count)
      : is_finished_(false), is_processing_count_(0) {
    for (size_t i = 0; i < thread_count; ++i) {
      workers_.emplace_back(std::async(std::launch::async, [&] {
        while (true) {
          std::future<F> future_to_process;
          {
            std::unique_lock<std::mutex> lock(futures_mutex_);
            futures_condition_.wait(
                lock, [&] { return is_finished_ || !futures_.empty(); });
            if (is_finished_ && futures_.empty()) { // nothing left to do, exit
              return;
            }
            future_to_process = std::move(*futures_.begin());
            futures_.pop_front();
            is_processing_count_++;
          }
          future_to_process.wait();
          {
            std::lock_guard<std::mutex> lock(futures_mutex_);
            is_processing_count_--;
          }
          processing_condition_.notify_one();
          {
            std::lock_guard<std::mutex> lock(done_futures_mutex_);
            done_futures_.emplace_back(std::move(future_to_process));
          }
        }
      }));
    }
  }

  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(futures_mutex_);
      if (is_finished_) {
        return;
      }
    }
    get();
  }

  template <typename T> void push(T lambda) {
    push(std::async(std::launch::async, lambda));
  }

  std::list<std::future<F>> &&get() {
    {
      std::lock_guard<std::mutex> lock(futures_mutex_);
      is_finished_ = true;
    }
    futures_condition_.notify_all();
    for (auto &worker : workers_) {
      worker.get(); // get any exceptions
    }
    return std::move(done_futures_);
  }

  void wait() {
    std::unique_lock<std::mutex> lock(futures_mutex_);
    processing_condition_.wait(
        lock, [&] { return is_processing_count_ == 0 && futures_.empty(); });
  }

private:
  void push(std::future<F> &&future) {
    {
      std::lock_guard<std::mutex> lock(futures_mutex_);
      if (is_finished_) {
        throw std::runtime_error(
            "Cannot push onto a thread pool after get() is used.");
      }
      futures_.emplace_back(std::move(future));
    }
    futures_condition_.notify_one();
  }

  std::vector<std::future<void>> workers_;
  std::list<std::future<F>> futures_;
  std::list<std::future<F>> done_futures_;
  std::mutex futures_mutex_;
  std::mutex done_futures_mutex_;
  std::condition_variable futures_condition_;
  std::condition_variable processing_condition_;
  bool is_finished_;
  size_t is_processing_count_;
};

} // namespace mcpd3
