# Progress Log

## 2026-07-13 22:00 PDT - Compact capacities with widened solver state

- Reproduced the accepted oracle-perturbed quantum-continuation experiments:
  128x128 reached objective `38393` in 2594 DD rounds and about 1.66 s;
  256x256 reached objective `155182` in 3973 DD rounds and about 3.3-3.6 s.
- Bisected the first rejected revision to checked-arithmetic commit `44bbbac`.
  UBSan on its parent found signed overflow while accumulating two valid arc
  flows into one 32-bit node balance: `-249043770 + -1960811640`.
- Kept source capacities, arc flows, and BK arc residuals in compact
  `Capacity`, while widening node balances, terminal residuals, Lagrange
  multipliers, and their persistent state to `Objective`.
- Added regressions for aggregate node balance beyond `Capacity`, alpha beyond
  `Capacity`, mixed-width BK augmentation, memory accounting, and mcpd4
  stateless/delta transport of widened alpha values.
- Bumped the streaming warm-state format to version 3 and mcpd4's transport
  protocol to version 9.
- Correct arithmetic changes the deterministic DD trajectory. The corrected
  Release runs reached the same optima at 3140 DD rounds / 5.51 s (128x128)
  and 4582 DD rounds / 22.40 s (256x256). The 128x128 corrected run also
  completed under UBSan with no signed-overflow report.
- Passed mcpd3's two tests in 32/64/128/GMP modes, mcpd4's five-test 32-bit
  integration suite plus widened transport tests in every precision mode, and
  all seven phase-unwrapping Release integration tests.

## 2026-07-13 15:47 PDT - Ratio-aware persistent flow refresh

- Extended primal-dual and dual-decomposition capacity refresh APIs with an
  optional widened rational flow scale.
- Require every internal forward/reverse capacity to have the supplied ratio
  before scaling; update signed arc flow with checked truncation toward zero,
  recompute node balances, and rebuild BK residual state.
- Validate all local partitions before applying a scaled DD refresh, so a
  rejection leaves every partition unchanged.
- Preserve alpha, last-alpha, and momentum independently from local flow.
- Added tests for both flow signs, fractional ratios, reset precedence,
  non-positive ratios, arithmetic overflow, non-proportional capacities,
  partition forwarding, balance reconstruction, and unchanged dual state.
- Passed both mcpd3 tests in all four capacity modes: 8/8 tests.

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

## 2026-07-13 15:08 PDT - Capacity precision performance validation

- Benchmarked pre-refactor 32-bit commit `43e0ada` against optimized precision
  commit `e6030dd` in isolated Release builds.
- Removed per-capacity Boost parsing from native DIMACS input and text-based
  integer conversion from the Lagrange hot path while retaining range checks.
- Current 32-bit total wall changed by -2.0% to +0.3% across PU 64x64,
  Waterloo bunny, and Waterloo gargoyle; no material regression was observed.
- Measured and documented 64-bit, 128-bit, and GMP runtime and peak-RSS costs
  in `CAPACITY_PERFORMANCE.md`.
- Verified identical objectives and iteration counts in all modes. Bunny
  matches its published objective directly; gargoyle's normalized objective
  plus its reported terminal-imbalance constant matches its published value.
- Passed both CTest tests in all four precision modes: 8/8 mode/test
  combinations.

## 2026-07-13 23:23 PDT - Deterministic historical DD replay

- Added opt-in `MCPD_LEGACY_32BIT_DD_REPLAY` for benchmark compatibility. It
  is restricted to 32-bit capacity builds and leaves the checked/widened
  production policy unchanged by default.
- Recreated the historical narrow Lagrange, node-balance, and BK terminal
  residual domains with explicit unsigned-bit wrapping. The replay therefore
  has deterministic two's-complement behavior without invoking signed-overflow
  undefined behavior.
- Added branch tests for replay type widths, maximum-to-minimum and
  minimum-to-maximum wrapping, wrapped aggregate node balance, normal checked
  overflow, memory estimates, and warm-state serialization types.
- Both normal and replay test builds pass. The replay tests also pass under
  UBSan with `-fno-sanitize-recover=undefined`; 64-bit replay configuration is
  rejected at configure time.
- Current phase code linked against this replay policy exactly reproduces the
  historical per-cut DD trajectories at 128x128 and 256x256. This is a
  diagnostic compatibility mode, not a production arithmetic policy.
