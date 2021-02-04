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

#include <unordered_set>
#include <map>
#include <vector>

namespace mcpd3 {

  class CycleCountingList {
    public:
      void addNode(long node) {
        list_.emplace_back(node);
        auto [iter,success] = seen_.insert(node);
        if ( !success ) { // seen before
          long current = node;
          size_t i = list_.size() - 2;
          std::vector<long> cycle;
          cycle.push_back(node);
          do {
            current = list_[i--];
            cycle.push_back(current);
          } while ( current != node );
          cycle_histogram_[cycle]++;
        }
      }

      int getMaxCycleCount() const {
        int max_count = 0;
        for (const auto &[cycle,count] : cycle_histogram_ ) {
          max_count = std::max(count,max_count);
        }
        return max_count;
      }

    private:
      std::vector<long> list_;
      std::unordered_set<long> seen_;
      std::map<std::vector<long>,int> cycle_histogram_;
  };

}
