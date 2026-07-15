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

## 2026-07-14 11:57 PDT - Cumulative scaled-epsilon regularization

- Replaced production binary sink anchors with persistent per-local-boundary
  weights. Every changed-alpha solve whose previous local label is sink adds
  the current epsilon to that weight; weights never decrease during a solve.
- Total regularization budget is now the sum of all persistent weights. Thus a
  continuing binary boundary disagreement consumes at least one epsilon unit
  per active regularized round and must eventually agree or reach the strict
  budget limit, which triggers objective-scale promotion or an over-budget
  result.
- Objective promotion now implements `factor * F + R`: primary capacities,
  alphas, and explicit arc flow scale, cumulative regularization remains
  unscaled, and the BK residual is reconstructed from the retained warm flow.
- Added the exact agreement certificate: when local labels agree and total
  regularization is strictly below one primary quantum, the feasible primary
  objective is itself the certified lower bound. During disagreement the
  conservative `regularized objective - total budget` certificate remains.
- Streaming warm state version 4 persists full cumulative weights and reads
  version 3 binary anchors compatibly. Tests cover repeated accumulation,
  source-label persistence, cycle-by-cycle budget growth, eviction/reload,
  strict promotion, and agreement certification.
- Passed complete CTest suites in 32-, 64-, 128-bit, GMP, and opt-in legacy
  replay builds. Legacy replay intentionally retains its historical binary
  anchor semantics.

## 2026-07-14 13:12 PDT - Disagreement-plateau regularization

- Added opt-in `DISAGREEMENT_PLATEAU_EPSILON` scheduling for native DD. Each
  scale starts with no new epsilon accumulation, waits for a separately
  configured number of iterations without a lower disagreement count, then
  enables cumulative regularization at fixed strength 1.
- Activating regularization resets lower-bound patience. If agreement is not
  reached during that normal patience window, the solver advances to the next
  scale and starts a fresh disagreement window. Existing cumulative weights
  remain part of the warm solver state and strict-budget promotion remains the
  exactness safeguard.
- Added branch tests for invalid observations/options, full plateau windows,
  disagreement improvement, scale reset, one-time activation, fixed unit
  strength, and a real two-partition plateau that consumes regularization
  budget while remaining strictly below the objective quantum.
- Passed complete CTest suites in 32-, 64-, 128-bit, GMP, and legacy replay
  builds.

## 2026-07-15 15:08 PDT - Opt-in edge-flow heat and weighted METIS

- Added opt-in per-original-edge counters to `PrimalDualMinCutSolver`. A count
  increments when a completed local maxflow changes the edge's stored net flow
  by a nonzero amount. Disabled tracking allocates no per-edge storage.
- Added global reconstruction through `DualDecomposition::arc_locations_` and
  retained partition labels only for tracking or explicitly weighted runs.
  Each crossing edge is read from its sole owning partition; cloned boundary
  nodes and terminal capacities cannot duplicate its count.
- Added validated positive METIS edge weights. Explicit weights are rejected
  by non-METIS partitioners and builds instead of being silently ignored.
- Tests cover tracking disabled/enabled/reset behavior, a reversed crossing
  edge represented exactly once, weighted partition selection, wrong-sized,
  zero, and out-of-range weights, and non-METIS rejection.
- Passed both Release 64-bit suites: 2/2 with METIS and 2/2 without METIS.

## 2026-07-15 - mcpd3-nh halo worktree baseline

- Created branch `exp/mcpd3-nh-halo` in the isolated worktree
  `/home/matt/software/experiments/mcpd3-nh-halo` from productized mcpd3
  commit `00ef729`.
- Wrote `HALO_DECOMPOSITION.md` to pin h1 compatibility, finite BFS halo
  membership, node-only consensus, exact edge/node multiplicity scaling, and
  verification gates before implementation.
- Configured a clean Release test build with the machine's validated local
  Boost headers. The unmodified baseline passed 2/2 CTest tests.
- Added a pure halo-layout planner under test-driven development. The tests
  first failed because the planner header did not exist, then passed after the
  implementation.
- The planner preserves historical one-owner h1 placement, computes finite
  multi-source BFS memberships for h2+, supports an explicit infinite-halo
  sentinel, intersects endpoint memberships for edge copies, and computes a
  checked LCM across every node and edge multiplicity.
- Tests cover exact h1/h2/infinite memberships, mixed two-way/three-way
  overlap requiring multiplier six, invalid depth/labels/endpoints, and LCM
  overflow. The full Release suite remains 2/2.
