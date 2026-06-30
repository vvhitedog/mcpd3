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

#include <cassert>
#include <chrono>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <graph/mcgraph.h>

namespace mcpd3 {

namespace _dimacs_implementation {

inline bool progress_enabled() {
  const char *value = std::getenv("MCPD3_PROGRESS");
  return value != nullptr && value[0] != '\0' && std::string(value) != "0";
}

inline void progress_report(const char *stage, long done, long total,
                            std::chrono::steady_clock::time_point start) {
  if (!progress_enabled()) {
    return;
  }
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
  const double pct = total > 0 ? 100.0 * static_cast<double>(done) / total : 0;
  const double rate = elapsed > 0 ? static_cast<double>(done) / elapsed : 0;
  const double eta = rate > 0 && total > done
                         ? static_cast<double>(total - done) / rate
                         : 0;
  std::fprintf(stderr,
               "mcpd3_progress stage=%s done=%ld total=%ld pct=%.2f "
               "elapsed_sec=%.1f eta_sec=%.1f rate=%.0f_per_sec\n",
               stage, done, total, pct, elapsed, eta, rate);
  std::fflush(stderr);
}

inline const char *skip_space(const char *p) {
  while (*p != '\0' && std::isspace(static_cast<unsigned char>(*p))) {
    ++p;
  }
  return p;
}

inline const char *skip_token(const char *p) {
  p = skip_space(p);
  while (*p != '\0' && !std::isspace(static_cast<unsigned char>(*p))) {
    ++p;
  }
  return p;
}

inline bool parse_int_token(const char *&p, int &value) {
  p = skip_space(p);
  bool negative = false;
  if (*p == '-') {
    negative = true;
    ++p;
  }
  if (!std::isdigit(static_cast<unsigned char>(*p))) {
    return false;
  }
  int parsed = 0;
  while (std::isdigit(static_cast<unsigned char>(*p))) {
    parsed = parsed * 10 + (*p - '0');
    ++p;
  }
  value = negative ? -parsed : parsed;
  return true;
}

inline bool parse_char_token(const char *&p, char &value) {
  p = skip_space(p);
  if (*p == '\0') {
    return false;
  }
  value = *p++;
  return true;
}

int remap_index(int original_index, int source_index, int sink_index) {
  int new_index = original_index;
  if (original_index > source_index) {
    --new_index;
  }
  if (original_index > sink_index) {
    --new_index;
  }
  return --new_index; // indices in DIMACS start at 1 but we wish for them
                      // to start at 0
}

template <typename ArcOperator, typename TerminalOperator>
void read_dimacs_general(const std::string &filename, ArcOperator arc_op,
                         TerminalOperator term_op) {
  const int line_length = 1024;
  char line[line_length];
  int n, m, s, t, cap, source, sink, _s, _t;
  char c;
  FILE *stream = nullptr;
  source = -1;
  sink = -1;
  if (!(stream = fopen(filename.c_str(), "r"))) {
    std::string err_msg = "failed to open file for reading: " + filename;
    throw std::runtime_error(err_msg.c_str());
  }
  while (fgets(line, line_length, stream)) {
    switch (*line) {
    case 'p':
      {
      const char *p = line + 1;
      p = skip_token(p); // max
      if (!parse_int_token(p, n) || !parse_int_token(p, m)) {
        std::string err_msg =
            "p line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      }
      break;
    case 'a':
      if (source == -1 || sink == -1) {
        std::string err_msg =
            "'a' line occured beforce setting source/sink in DIMACS file:" +
            filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      {
      const char *p = line + 1;
      if (!parse_int_token(p, s) || !parse_int_token(p, t) ||
          !parse_int_token(p, cap)) {
        std::string err_msg =
            "'a' line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      }
      if (t == source || s == sink) {
        std::string err_msg =
            "specified source or sink as target or source node incorrectly:" +
            filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      if (s != source && t != sink) { // handle non terminal arc
        _s = remap_index(s, source, sink);
        _t = remap_index(t, source, sink);
        arc_op(_s, _t, cap);
      } else { // handle terminal arc
        if (s == source) {
          _t = remap_index(t, source, sink);
          term_op(true, _t, cap);
        } else if (t == sink) {
          _s = remap_index(s, source, sink);
          term_op(false, _s, cap);
        }
      }
      break;
    case 'n':
      {
      const char *p = line + 1;
      if (!parse_int_token(p, s) || !parse_char_token(p, c)) {
        std::string err_msg =
            "'n' line is malformed in DIMACS file:" + filename + "\n";
        throw std::runtime_error(err_msg.c_str());
      }
      }
      if (c == 's') {
        source = s;
      }
      if (c == 't') {
        sink = s;
      }
      break;
    case 'c':
      // comment
      break;
    default:
      break;
    }
  }
  fclose(stream);
}
} // namespace _dimacs_implementation

MinCutGraph read_dimacs(const std::string &filename) {
  MinCutGraph g;
  g.nnode = 0;
  std::unordered_map<int, std::unordered_map<int, int>> arc_adjacency;

  auto arc_op = [&](int s, int t, int cap) {
    g.nnode = std::max(g.nnode, s + 1);
    g.nnode = std::max(g.nnode, t + 1);
    if (arc_adjacency[s].find(t) == arc_adjacency[s].end()) {
      arc_adjacency[s][t] = cap;
      if (arc_adjacency[t].find(s) == arc_adjacency[t].end()) {
        arc_adjacency[t][s] = 0;
      }
    } else {
      arc_adjacency[s][t] += cap;
    }
  };

  long imbalance = 0;
  auto term_op = [&](bool is_source, int n, int cap) {
    g.nnode = std::max(g.nnode, n + 1);
    g.terminal_capacities.resize(g.nnode, 0);
    if (is_source) {
      auto old_cap = g.terminal_capacities[n];
      if (old_cap < 0) {
        imbalance += std::min(-old_cap, cap);
      }
      g.terminal_capacities[n] += cap;
    } else {
      auto old_cap = g.terminal_capacities[n];
      if (old_cap > 0) {
        imbalance += std::min(old_cap, cap);
      }
      g.terminal_capacities[n] -= cap;
    }
  };

  _dimacs_implementation::read_dimacs_general(filename, arc_op, term_op);

  if (imbalance > 0) {
    printf("WARNING: imbalance when reading dimacs graph: %lu\n", imbalance);
  }

  // process arcs after the fact
  g.narc = 0;
  for (const auto &[s, arc_list] : arc_adjacency) {
    for (const auto &[t, cap] : arc_list) {
      if (s < t) {
        g.arc_capacities.push_back(cap);
        g.arc_capacities.push_back(arc_adjacency[t][s]);
        g.arcs.push_back(s);
        g.arcs.push_back(t);
        g.narc++;
      }
    }
  }
  g.terminal_capacities.resize(g.nnode, 0); // ensure size
  g.arcs.shrink_to_fit();
  g.terminal_capacities.shrink_to_fit();
  g.arc_capacities.shrink_to_fit();
  return std::move(g);
}

MinCutGraph read_dimacs_directed_streaming(const std::string &filename) {
  MinCutGraph g;
  g.nnode = 0;
  g.narc = 0;

  auto arc_op = [&](int s, int t, int cap) {
    g.nnode = std::max(g.nnode, s + 1);
    g.nnode = std::max(g.nnode, t + 1);
    g.arc_capacities.push_back(cap);
    g.arc_capacities.push_back(0);
    g.arcs.push_back(s);
    g.arcs.push_back(t);
    ++g.narc;
  };

  long imbalance = 0;
  auto term_op = [&](bool is_source, int n, int cap) {
    g.nnode = std::max(g.nnode, n + 1);
    g.terminal_capacities.resize(g.nnode, 0);
    if (is_source) {
      auto old_cap = g.terminal_capacities[n];
      if (old_cap < 0) {
        imbalance += std::min(-old_cap, cap);
      }
      g.terminal_capacities[n] += cap;
    } else {
      auto old_cap = g.terminal_capacities[n];
      if (old_cap > 0) {
        imbalance += std::min(old_cap, cap);
      }
      g.terminal_capacities[n] -= cap;
    }
  };

  _dimacs_implementation::read_dimacs_general(filename, arc_op, term_op);

  if (imbalance > 0) {
    printf("WARNING: imbalance when reading dimacs graph: %lu\n", imbalance);
  }

  g.terminal_capacities.resize(g.nnode, 0);
  g.arcs.shrink_to_fit();
  g.terminal_capacities.shrink_to_fit();
  g.arc_capacities.shrink_to_fit();
  return std::move(g);
}

MinCutGraph read_dimacs_symmetric_streaming(const std::string &filename) {
  const int line_length = 1024;
  char line[line_length];
  int declared_nodes = 0;
  int declared_arcs = 0;
  int source = -1;
  int sink = -1;

  {
    FILE *stream = fopen(filename.c_str(), "r");
    if (!stream) {
      throw std::runtime_error("failed to open file for reading: " + filename);
    }
    while (fgets(line, line_length, stream)) {
      switch (*line) {
      case 'p': {
        const char *p = line + 1;
        p = _dimacs_implementation::skip_token(p); // max
        if (!_dimacs_implementation::parse_int_token(p, declared_nodes) ||
            !_dimacs_implementation::parse_int_token(p, declared_arcs)) {
          fclose(stream);
          throw std::runtime_error("malformed p line in DIMACS file: " +
                                   filename);
        }
        break;
      }
      case 'n': {
        int node = 0;
        char terminal = '\0';
        const char *p = line + 1;
        if (!_dimacs_implementation::parse_int_token(p, node) ||
            !_dimacs_implementation::parse_char_token(p, terminal)) {
          fclose(stream);
          throw std::runtime_error("malformed n line in DIMACS file: " +
                                   filename);
        }
        if (terminal == 's') {
          source = node;
        } else if (terminal == 't') {
          sink = node;
        }
        break;
      }
      default:
        break;
      }
      if (declared_nodes > 0 && declared_arcs > 0 && source > 0 && sink > 0) {
        break;
      }
    }
    fclose(stream);
  }

  if (declared_nodes <= 0 || declared_arcs <= 0 || source <= 0 || sink <= 0) {
    throw std::runtime_error("missing DIMACS header/source/sink in file: " +
                             filename);
  }

  long internal_pair_count = 0;
  long terminal_arc_count = 0;
  bool has_pending = false;
  int pending_s = 0;
  int pending_t = 0;
  int pending_cap = 0;
  const bool report_progress = _dimacs_implementation::progress_enabled();
  const long progress_interval = 10000000;

  {
    auto pass_start = std::chrono::steady_clock::now();
    long arcs_seen = 0;
    FILE *stream = fopen(filename.c_str(), "r");
    if (!stream) {
      throw std::runtime_error("failed to open file for reading: " + filename);
    }
    while (fgets(line, line_length, stream)) {
      if (*line != 'a') {
        continue;
      }
      ++arcs_seen;
      if (report_progress && arcs_seen % progress_interval == 0) {
        _dimacs_implementation::progress_report(
            "dimacs_stream_count", arcs_seen, declared_arcs, pass_start);
      }
      int s = 0;
      int t = 0;
      int cap = 0;
      const char *p = line + 1;
      if (!_dimacs_implementation::parse_int_token(p, s) ||
          !_dimacs_implementation::parse_int_token(p, t) ||
          !_dimacs_implementation::parse_int_token(p, cap)) {
        fclose(stream);
        throw std::runtime_error("malformed a line in DIMACS file: " +
                                 filename);
      }
      if (t == source || s == sink) {
        fclose(stream);
        throw std::runtime_error(
            "DIMACS source/sink appears in invalid arc orientation: " +
            filename);
      }
      const bool terminal_arc = (s == source || t == sink);
      if (terminal_arc) {
        if (has_pending) {
          fclose(stream);
          throw std::runtime_error(
              "terminal arc encountered between an internal arc and its "
              "adjacent reverse arc in DIMACS file: " +
              filename);
        }
        ++terminal_arc_count;
        continue;
      }
      if (!has_pending) {
        pending_s = s;
        pending_t = t;
        pending_cap = cap;
        has_pending = true;
      } else {
        if (s != pending_t || t != pending_s) {
          fclose(stream);
          throw std::runtime_error(
              "non-terminal arcs are not adjacent reverse pairs in DIMACS "
              "file: " +
              filename);
        }
        (void)pending_cap;
        ++internal_pair_count;
        has_pending = false;
      }
    }
    fclose(stream);
    _dimacs_implementation::progress_report("dimacs_stream_count", arcs_seen,
                                            declared_arcs, pass_start);
  }
  if (has_pending) {
    throw std::runtime_error("DIMACS file ended with an unmatched internal arc: " +
                             filename);
  }

  MinCutGraph g;
  g.nnode = declared_nodes - 2;
  g.narc = static_cast<int>(internal_pair_count);
  g.arcs.reserve(static_cast<size_t>(2) * internal_pair_count);
  g.arc_capacities.reserve(static_cast<size_t>(2) * internal_pair_count);
  g.terminal_capacities.resize(g.nnode, 0);

  long imbalance = 0;
  has_pending = false;
  auto build_start = std::chrono::steady_clock::now();
  long arcs_seen = 0;
  FILE *stream = fopen(filename.c_str(), "r");
  if (!stream) {
    throw std::runtime_error("failed to open file for reading: " + filename);
  }
  while (fgets(line, line_length, stream)) {
    if (*line != 'a') {
      continue;
    }
    ++arcs_seen;
    if (report_progress && arcs_seen % progress_interval == 0) {
      _dimacs_implementation::progress_report(
          "dimacs_stream_build", arcs_seen, declared_arcs, build_start);
    }
    int s = 0;
    int t = 0;
    int cap = 0;
    const char *p = line + 1;
    if (!_dimacs_implementation::parse_int_token(p, s) ||
        !_dimacs_implementation::parse_int_token(p, t) ||
        !_dimacs_implementation::parse_int_token(p, cap)) {
      fclose(stream);
      throw std::runtime_error("malformed a line in DIMACS file: " + filename);
    }
    const bool terminal_arc = (s == source || t == sink);
    if (terminal_arc) {
      const int n = (s == source)
                        ? _dimacs_implementation::remap_index(t, source, sink)
                        : _dimacs_implementation::remap_index(s, source, sink);
      if (s == source) {
        auto old_cap = g.terminal_capacities[n];
        if (old_cap < 0) {
          imbalance += std::min(-old_cap, cap);
        }
        g.terminal_capacities[n] += cap;
      } else if (t == sink) {
        auto old_cap = g.terminal_capacities[n];
        if (old_cap > 0) {
          imbalance += std::min(old_cap, cap);
        }
        g.terminal_capacities[n] -= cap;
      }
      continue;
    }

    if (!has_pending) {
      pending_s = s;
      pending_t = t;
      pending_cap = cap;
      has_pending = true;
    } else {
      const int remapped_s =
          _dimacs_implementation::remap_index(pending_s, source, sink);
      const int remapped_t =
          _dimacs_implementation::remap_index(pending_t, source, sink);
      g.arcs.push_back(remapped_s);
      g.arcs.push_back(remapped_t);
      g.arc_capacities.push_back(pending_cap);
      g.arc_capacities.push_back(cap);
      has_pending = false;
    }
  }
  fclose(stream);
  _dimacs_implementation::progress_report("dimacs_stream_build", arcs_seen,
                                          declared_arcs, build_start);

  if (imbalance > 0) {
    printf("WARNING: imbalance when reading dimacs graph: %lu\n", imbalance);
  }
  return std::move(g);
}

} // namespace mcpd3
