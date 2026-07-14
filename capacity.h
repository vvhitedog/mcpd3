#pragma once

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

#if defined(MCPD_CAPACITY_MODE_GMP)
#include <gmpxx.h>
#endif

#if (defined(MCPD_CAPACITY_MODE_32) + defined(MCPD_CAPACITY_MODE_64) +         \
     defined(MCPD_CAPACITY_MODE_128) + defined(MCPD_CAPACITY_MODE_GMP)) > 1
#error "Only one MCPD capacity mode may be selected"
#endif

#if !defined(MCPD_CAPACITY_MODE_32) && !defined(MCPD_CAPACITY_MODE_64) &&     \
    !defined(MCPD_CAPACITY_MODE_128) && !defined(MCPD_CAPACITY_MODE_GMP)
#define MCPD_CAPACITY_MODE_32 1
#endif

namespace mcpd3 {

#if defined(MCPD_CAPACITY_MODE_32)
using Capacity = std::int32_t;
using Objective = std::int64_t;
#elif defined(MCPD_CAPACITY_MODE_64)
using Capacity = std::int64_t;
using Objective = __int128_t;
#elif defined(MCPD_CAPACITY_MODE_128)
using Capacity = __int128_t;
using Objective = boost::multiprecision::int256_t;
#else
using Capacity = mpz_class;
using Objective = mpz_class;
#endif

// Source capacities stay compact; solver potentials must hold their sums.
using Lagrange = Objective;

inline const char *capacity_mode_name() {
#if defined(MCPD_CAPACITY_MODE_32)
  return "32";
#elif defined(MCPD_CAPACITY_MODE_64)
  return "64";
#elif defined(MCPD_CAPACITY_MODE_128)
  return "128";
#else
  return "gmp";
#endif
}

inline constexpr unsigned capacity_storage_bits() {
#if defined(MCPD_CAPACITY_MODE_32)
  return 32;
#elif defined(MCPD_CAPACITY_MODE_64)
  return 64;
#elif defined(MCPD_CAPACITY_MODE_128)
  return 128;
#else
  return 0;
#endif
}

inline constexpr bool capacity_is_bounded() {
  return capacity_storage_bits() != 0;
}

inline std::string int128_to_string(__int128_t value) {
  if (value == 0) {
    return "0";
  }
  const bool negative = value < 0;
  __uint128_t magnitude = negative
                                ? static_cast<__uint128_t>(-(value + 1)) + 1
                                : static_cast<__uint128_t>(value);
  std::string result;
  while (magnitude != 0) {
    result.push_back(static_cast<char>('0' + magnitude % 10));
    magnitude /= 10;
  }
  if (negative) {
    result.push_back('-');
  }
  std::reverse(result.begin(), result.end());
  return result;
}

template <typename Integer>
inline std::string integer_to_string(const Integer &value) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if constexpr (std::is_constructible_v<mpz_class, Integer>) {
    return mpz_class(value).get_str();
  } else
#endif
      if constexpr (std::is_same_v<std::decay_t<Integer>, __int128_t>) {
    return int128_to_string(value);
  } else if constexpr (boost::multiprecision::is_number<Integer>::value) {
    return value.template convert_to<std::string>();
  } else {
    return std::to_string(value);
  }
}

template <typename Integer>
inline double integer_to_double(const Integer &value) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if constexpr (std::is_constructible_v<mpz_class, Integer>) {
    return mpz_class(value).get_d();
  } else
#endif
      if constexpr (boost::multiprecision::is_number<Integer>::value) {
    return value.template convert_to<double>();
  } else {
    return static_cast<double>(value);
  }
}

template <typename Integer>
inline Integer parse_bounded_integer(const std::string &text) {
  if (text.empty()) {
    throw std::invalid_argument("integer text is empty");
  }
  boost::multiprecision::cpp_int parsed;
  try {
    parsed = boost::multiprecision::cpp_int(text);
  } catch (const std::exception &) {
    throw std::invalid_argument("invalid decimal integer: " + text);
  }
  const boost::multiprecision::cpp_int minimum(
      integer_to_string(std::numeric_limits<Integer>::min()));
  const boost::multiprecision::cpp_int maximum(
      integer_to_string(std::numeric_limits<Integer>::max()));
  if (parsed < minimum || parsed > maximum) {
    throw std::out_of_range("decimal integer is outside configured range");
  }
  return parsed.template convert_to<Integer>();
}

template <typename Integer>
inline bool parse_bounded_integer_chars(const char *begin, const char *end,
                                        Integer &value) {
  if (begin == end) {
    return false;
  }

  bool negative = false;
  if (*begin == '-' || *begin == '+') {
    negative = *begin == '-';
    ++begin;
  }
  if (begin == end) {
    return false;
  }

  using Unsigned = std::make_unsigned_t<Integer>;
  const Unsigned positive_limit =
      static_cast<Unsigned>(std::numeric_limits<Integer>::max());
  const Unsigned limit =
      negative ? static_cast<Unsigned>(positive_limit + 1) : positive_limit;
  Unsigned magnitude = 0;
  for (const char *cursor = begin; cursor != end; ++cursor) {
    if (*cursor < '0' || *cursor > '9') {
      return false;
    }
    const Unsigned digit = static_cast<Unsigned>(*cursor - '0');
    if (digit > limit || magnitude > (limit - digit) / 10) {
      return false;
    }
    magnitude = static_cast<Unsigned>(magnitude * 10 + digit);
  }

  if (negative) {
    if (magnitude == static_cast<Unsigned>(positive_limit + 1)) {
      value = std::numeric_limits<Integer>::min();
    } else {
      value = static_cast<Integer>(-static_cast<Integer>(magnitude));
    }
  } else {
    value = static_cast<Integer>(magnitude);
  }
  return true;
}

inline bool parse_capacity_chars(const char *begin, const char *end,
                                 Capacity &value) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if (begin == end) {
    return false;
  }
  const char *cursor = begin;
  const char *parse_begin = begin;
  if (*cursor == '-' || *cursor == '+') {
    if (*cursor == '+') {
      parse_begin = cursor + 1;
    }
    ++cursor;
  }
  if (cursor == end) {
    return false;
  }
  for (; cursor != end; ++cursor) {
    if (*cursor < '0' || *cursor > '9') {
      return false;
    }
  }
  try {
    value = Capacity(std::string(parse_begin, end));
    return true;
  } catch (const std::exception &) {
    return false;
  }
#elif defined(MCPD_CAPACITY_MODE_32) || defined(MCPD_CAPACITY_MODE_64)
  if (begin == end) {
    return false;
  }
  if (*begin == '+') {
    ++begin;
  }
  if (begin == end) {
    return false;
  }
  const auto result = std::from_chars(begin, end, value, 10);
  return result.ec == std::errc{} && result.ptr == end;
#else
  return parse_bounded_integer_chars(begin, end, value);
#endif
}

inline Capacity parse_capacity(const std::string &text) {
  Capacity value = 0;
  if (!parse_capacity_chars(text.data(), text.data() + text.size(), value)) {
    throw std::invalid_argument("invalid or out-of-range decimal capacity: " +
                                text);
  }
  return value;
}

template <typename Integer,
          std::enable_if_t<std::is_integral_v<std::decay_t<Integer>>, int> = 0>
inline Capacity capacity_from_integer(Integer value) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if constexpr (std::is_signed_v<Integer> && sizeof(Integer) <= sizeof(long)) {
    return Capacity(static_cast<long>(value));
  } else if constexpr (std::is_unsigned_v<Integer> &&
                       sizeof(Integer) <= sizeof(unsigned long)) {
    return Capacity(static_cast<unsigned long>(value));
  } else {
    return parse_capacity(integer_to_string(value));
  }
#else
  if constexpr (std::is_signed_v<Integer>) {
    if constexpr (std::numeric_limits<Integer>::digits <=
                  std::numeric_limits<Capacity>::digits) {
      return static_cast<Capacity>(value);
    } else {
      if (value < static_cast<Integer>(std::numeric_limits<Capacity>::min()) ||
          value > static_cast<Integer>(std::numeric_limits<Capacity>::max())) {
        throw std::overflow_error("integer does not fit capacity type");
      }
      return static_cast<Capacity>(value);
    }
  } else {
    if constexpr (std::numeric_limits<Integer>::digits <=
                  std::numeric_limits<Capacity>::digits) {
      return static_cast<Capacity>(value);
    } else {
      if (value > static_cast<Integer>(std::numeric_limits<Capacity>::max())) {
        throw std::overflow_error("integer does not fit capacity type");
      }
      return static_cast<Capacity>(value);
    }
  }
#endif
}

template <typename Integer,
          std::enable_if_t<!std::is_integral_v<std::decay_t<Integer>>, int> = 0>
inline Capacity capacity_from_integer(const Integer &value) {
  return parse_capacity(integer_to_string(value));
}

template <typename Integer>
inline std::vector<Capacity>
capacity_vector_from(const std::vector<Integer> &values) {
  std::vector<Capacity> result;
  result.reserve(values.size());
  for (const auto &value : values) {
    result.push_back(capacity_from_integer(value));
  }
  return result;
}

inline Objective parse_objective(const std::string &text) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  try {
    return Objective(text);
  } catch (const std::exception &) {
    throw std::invalid_argument("invalid decimal objective: " + text);
  }
#else
  return parse_bounded_integer<Objective>(text);
#endif
}

template <typename Integer,
          std::enable_if_t<std::is_integral_v<std::decay_t<Integer>>, int> = 0>
inline Lagrange lagrange_from_integer(Integer value) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if constexpr (std::is_signed_v<Integer> && sizeof(Integer) <= sizeof(long)) {
    return Lagrange(static_cast<long>(value));
  } else if constexpr (std::is_unsigned_v<Integer> &&
                       sizeof(Integer) <= sizeof(unsigned long)) {
    return Lagrange(static_cast<unsigned long>(value));
  } else {
    return parse_objective(integer_to_string(value));
  }
#else
  if constexpr (std::is_signed_v<Integer>) {
    if constexpr (std::numeric_limits<Integer>::digits <=
                  std::numeric_limits<Lagrange>::digits) {
      return static_cast<Lagrange>(value);
    }
  } else if constexpr (std::numeric_limits<Integer>::digits <=
                       std::numeric_limits<Lagrange>::digits) {
    return static_cast<Lagrange>(value);
  }
  return parse_objective(integer_to_string(value));
#endif
}

template <typename Integer,
          std::enable_if_t<!std::is_integral_v<std::decay_t<Integer>>, int> = 0>
inline Lagrange lagrange_from_integer(const Integer &value) {
  return parse_objective(integer_to_string(value));
}

inline Objective widen_capacity(const Capacity &value) {
  return static_cast<Objective>(value);
}

inline Objective absolute_capacity(const Capacity &value) {
  const Objective widened = widen_capacity(value);
  return widened < 0 ? -widened : widened;
}

template <typename Integer>
inline constexpr bool integer_is_bounded() {
#if defined(MCPD_CAPACITY_MODE_GMP)
  if constexpr (std::is_same_v<std::decay_t<Integer>, mpz_class>) {
    return false;
  }
#endif
  return std::numeric_limits<Integer>::is_bounded;
}

template <typename Integer>
inline Integer checked_add(const Integer &lhs, const Integer &rhs,
                           const char *message = "integer addition overflow") {
  if constexpr (!integer_is_bounded<Integer>()) {
    return lhs + rhs;
  } else {
    if ((rhs > 0 && lhs > std::numeric_limits<Integer>::max() - rhs) ||
        (rhs < 0 && lhs < std::numeric_limits<Integer>::min() - rhs)) {
      throw std::overflow_error(message);
    }
    return lhs + rhs;
  }
}

template <typename Integer>
inline Integer checked_subtract(
    const Integer &lhs, const Integer &rhs,
    const char *message = "integer subtraction overflow") {
  if constexpr (!integer_is_bounded<Integer>()) {
    return lhs - rhs;
  } else {
    if ((rhs > 0 && lhs < std::numeric_limits<Integer>::min() + rhs) ||
        (rhs < 0 && lhs > std::numeric_limits<Integer>::max() + rhs)) {
      throw std::overflow_error(message);
    }
    return lhs - rhs;
  }
}

template <typename Integer>
inline Integer checked_scale(
    const Integer &value, long scale,
    const char *message = "integer multiplication overflow") {
  if (scale <= 0) {
    throw std::invalid_argument("scale factor must be positive");
  }
  if constexpr (!integer_is_bounded<Integer>()) {
    return value * scale;
  } else {
    if (value > 0 && value > std::numeric_limits<Integer>::max() / scale) {
      throw std::overflow_error(message);
    }
    if (value < 0 && value < std::numeric_limits<Integer>::min() / scale) {
      throw std::overflow_error(message);
    }
    return value * scale;
  }
}

template <typename Integer>
inline Integer checked_multiply_by_positive(
    const Integer &value, const Integer &multiplier,
    const char *message = "integer multiplication overflow") {
  if (multiplier <= 0) {
    throw std::invalid_argument("multiplier must be positive");
  }
  if constexpr (!integer_is_bounded<Integer>()) {
    return value * multiplier;
  } else {
    if (value > 0 &&
        value > std::numeric_limits<Integer>::max() / multiplier) {
      throw std::overflow_error(message);
    }
    if (value < 0 &&
        value < std::numeric_limits<Integer>::min() / multiplier) {
      throw std::overflow_error(message);
    }
    return value * multiplier;
  }
}

inline bool capacities_have_ratio(const Capacity &old_capacity,
                                  const Capacity &new_capacity,
                                  const Objective &numerator,
                                  const Objective &denominator) {
  if (numerator <= 0 || denominator <= 0) {
    throw std::invalid_argument("flow scale ratio must be positive");
  }
  return checked_multiply_by_positive(
             widen_capacity(old_capacity), numerator,
             "capacity ratio multiplication overflow") ==
         checked_multiply_by_positive(
             widen_capacity(new_capacity), denominator,
             "capacity ratio multiplication overflow");
}

inline Capacity checked_scale_capacity(const Capacity &value, long scale,
                                        bool saturate = false) {
  try {
    return checked_scale(value, scale, "capacity multiplication overflow");
  } catch (const std::overflow_error &) {
    if (!saturate || !capacity_is_bounded()) {
      throw;
    }
#if defined(MCPD_CAPACITY_MODE_GMP)
    throw;
#else
    return value < 0 ? std::numeric_limits<Capacity>::min()
                     : std::numeric_limits<Capacity>::max();
#endif
  }
}

inline Capacity narrow_objective_to_capacity(const Objective &value,
                                              bool saturate = false) {
#if defined(MCPD_CAPACITY_MODE_GMP)
  return value;
#else
  const Objective minimum =
      static_cast<Objective>(std::numeric_limits<Capacity>::min());
  const Objective maximum =
      static_cast<Objective>(std::numeric_limits<Capacity>::max());
  if (value < minimum || value > maximum) {
    if (saturate) {
      return value < 0 ? std::numeric_limits<Capacity>::min()
                       : std::numeric_limits<Capacity>::max();
    }
    throw std::overflow_error("objective does not fit capacity type");
  }
  return static_cast<Capacity>(value);
#endif
}

inline Capacity checked_scale_capacity_ratio(
    const Capacity &value, const Objective &numerator,
    const Objective &denominator) {
  if (numerator <= 0 || denominator <= 0) {
    throw std::invalid_argument("flow scale ratio must be positive");
  }
  const Objective product = checked_multiply_by_positive(
      widen_capacity(value), numerator, "flow scale multiplication overflow");
  return narrow_objective_to_capacity(product / denominator);
}

inline Capacity capacity_test_extreme_value() {
#if defined(MCPD_CAPACITY_MODE_GMP)
  Capacity value = 1;
  value <<= 521;
  return value + 123456789;
#else
  return std::numeric_limits<Capacity>::max();
#endif
}

} // namespace mcpd3
