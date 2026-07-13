# Progress Log

## 2026-07-13 11:37 PDT - Build-time capacity precision

- Added `MCPD_CAPACITY_MODE=32|64|128|gmp`; 32-bit remains the default.
- Added checked `Capacity` and widened `Objective` types throughout BK,
  primal-dual solving, dual decomposition, partition workers/coordinator,
  DIMACS input, local streaming state, examples, and the native benchmark.
- GMP is the arbitrary-precision backend. Nontrivial GMP graph storage uses
  constructed heap arrays; fixed-width modes retain mmap support.
- Added configured-extreme tests covering BK reallocation, decimal round trips,
  DIMACS, primal-dual solving, in-process packages, and disk-streamed packages.
- Built every mcpd3 target in all four modes.
- Passed both CTest tests in all four modes: 8/8 mode/test combinations.

## 2026-07-13 11:54 PDT - Shared precision build configuration

- Factored capacity-mode dependency and compile-definition setup into
  `cmake/McpdCapacity.cmake` so mcpd4 and PU can consume exactly the same
  precision contract as mcpd3 without duplicating it.

## 2026-07-13 12:55 PDT - Precision-safe legacy CSR storage

- Changed `CsrGraph` and `read_dimacs_to_csr` defaults from hardcoded
  `int/long` to the configured `Capacity/Objective` types.
- Generalized the mmap array API so trivially copyable fixed-width values stay
  file-backed while nontrivial GMP values receive constructed heap storage.
- Added a DIMACS-to-CSR extreme-value test; GMP carries a value above 521 bits
  through edge sorting, terminal accumulation, CSR storage, and cut evaluation.
- Made CSR capacity accumulation and cut evaluation checked, and corrected the
  partition-label array to use the node index type instead of the capacity
  type.

## 2026-07-13 13:00 PDT - Checked solver arithmetic audit

- Replaced unchecked signed arithmetic in primal-dual flow bookkeeping,
  residual construction, Lagrange aggregation, regularization accounting,
  warm-start reconstruction, and primal/mincut objective evaluation.
- Aggregate flow and cut totals now use the widened `Objective` domain while
  capacity-domain overflow raises a deterministic exception.
- Added a two-component extreme-capacity regression whose objective is twice
  the largest configured capacity value.
- Rebuilt every mcpd3 target and passed both CTest tests in all four precision
  modes: 8/8 mode/test combinations.
