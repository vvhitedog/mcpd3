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

## 2026-07-15 - Exact halo subproblem integration

- Added `DualDecompositionOptions::halo_depth`, defaulting to one, and exposed
  the checked halo objective multiplier. The h1 construction retains the old
  one-owner arc and terminal-location arrays; multi-location storage is only
  allocated for h2+.
- For h2+, local induced edges and copied unaries receive exact integer factors
  `Q_h / r_e` and `Q_h / r_v`. All duplicated nodes receive the existing signed
  pairwise consensus constraints. Original primal capacities, effective
  objective scale, and explicit regularization budget limits include `Q_h`.
- Generalized capacity replacement and flow-heat reconstruction to every
  local copy. Replacement reapplies the fixed multiplicity factors; objective
  promotion then scales the resulting local state through the existing path.
- Added a package objective multiplier so in-process/streaming worker
  coordinators normalize raw halo objectives correctly. Coordinators reject
  mismatched multipliers, and package validation rejects nonpositive values.
- Tests first failed on the absent halo options/package metadata. They now
  cover byte-equivalent h1 exports, exact h2 package topology/capacities,
  exhaustive objective equality on eight randomized directed graphs and all
  64 cuts at h1/h2/h3/infinite, direct objective certification, capacity
  replacement plus promotion, native/package coordinator equivalence, invalid
  metadata, and duplicate-copy flow heat/reset.
- The complete Release suite passes 2/2 after integration. `git diff --check`
  passes.

## 2026-07-15 16:23 PDT - Native mcpd3-nh selector

- Added `--halo-depth N|infinite` to both native dual-decomposition command
  lines. The default remains one, so existing mcpd3-n invocations retain the
  historical decomposition without an extra allocation or package copy.
- Both tools report the requested depth and the resulting exact halo objective
  multiplier. Invalid zero, negative, malformed, and overflowing finite depths
  are rejected rather than silently selecting another formulation.
- Documented the compatibility and experimental infinite modes in the README.
- Rebuilt the Release tree and passed the complete 2/2 CTest suite;
  `git diff --check` passes.
## 2026-07-15 17:10 PDT - Halo-depth-one performance equivalence

- Hoisted the depth-one/deeper-halo arc-distribution branch outside the
  original-edge loop. The depth-one path now executes the original compact
  one-owner construction loop without a per-edge halo condition.
- Compared `mcpd3-nh --halo-depth 1` against the pre-halo `mcpd3-n` binary on
  the physical 64x64 wrapped-ramp fixture using two affined physical CPU cores.
  Four interleaved 30-run blocks produced baseline medians of 32.66 and
  32.96 ms and halo-depth-one medians of 32.56 and 32.86 ms.
- Every run returned objective 63 with the identical seven-cut, 151-DD-round
  trajectory. Callgrind measured 202.76 million instructions for the halo
  binary versus 202.04 million for the baseline, a 0.36% difference.
## 2026-07-15 17:24 PDT - Deeper-halo regularization promotion coverage

- Added a deterministic h2 path fixture whose two local subproblems initially
  prefer opposite labels on four duplicated nodes.
- With primary objective scale 10 and halo multiplier `Q=2`, cumulative
  scaled-epsilon regularization reaches budget 40 against the strict effective
  limit 20. The solver promotes once, retains `Q`, reaches agreement at
  effective scale 200, and finishes below budget.
- The regression compares the final raw and normalized certificates against a
  direct whole-graph mincut, proving that promotion preserves exact halo
  scaling rather than only exercising the control-flow branch.
## 2026-07-15 18:16 PDT - Final halo validation matrix

- Rebuilt the final tests in Release mode and passed 2/2 CTest targets in
  capacity modes 32, 64, 128, and GMP, plus 32-bit historical replay and
  64-bit METIS configurations. This covers checked fixed-width arithmetic,
  arbitrary precision, the compatibility replay path, and both partitioners.
- The phase adapter and `mcpd3-nh` benchmark selector pass all 30 tests against
  this final core revision. `git diff --check` is clean in both worktrees.

## 2026-07-15 - Allocation-free incremental cut maintenance

- Added a forced-full-recompute oracle test that drives repeated whole-graph,
  checkerboard, grouped, and randomized alpha changes. It verifies incremental
  labels and objectives after every solve and explicitly exercises edges whose
  two endpoints change in the same round.
- Replaced the fresh per-solve edge `unordered_set` in incremental cut-value
  maintenance with one reusable byte per node. Changed labels remain immutable
  while edge deltas are evaluated, and the original source endpoint owns an
  edge when both endpoints changed.
- All Release tests pass in capacity modes 32, 64, 128, and GMP, historical
  32-bit replay, and 64-bit METIS. The phase integration passes 36/36 tests.
- Pinned three-run Ghiglia-Pritt profiles retained identical objectives, cut
  counts, and DD rounds. Spiral h1/h2/h5 improved from 6.325/4.484/10.766 s to
  3.954/3.200/7.985 s. Head h1/h5 improved from 12.062/11.143 s to
  7.550/8.472 s. The remaining prominent allocation cost comes from BK's
  separate changed-flow-arc hash, not cut-value maintenance.

## 2026-07-15 - Allocation-free changed-flow arc collection

- Preserved BK's existing changed-arc hash overload for compatibility and
  added a contiguous logical-edge collector for the primal-dual hot path.
  One byte per logical edge stores a generation mark, deduplicating sister
  arcs and repeated path updates without allocations. The generation wraps
  through a full clear once per 255 incremental solves.
- The full oracle test performs 360 solves, so it validates the generation
  wrap as well as randomized label and flow changes. Release tests pass in
  32, 64, 128, GMP, legacy replay, and METIS builds; phase tests pass 36/36.
- Final three-run wall times, relative to the original hash-based baseline,
  are Spiral h1 3.798 s (-40.0%), h2 2.906 s (-35.2%), and h5 7.074 s
  (-34.3%); Head h1 7.237 s (-40.0%) and h5 7.687 s (-31.0%). Objectives,
  cut counts, DD rounds, and promotions are unchanged.
- A post-change call profile has zero lost samples. Changed-arc hash insertion
  fell from 7.1% self time to zero; dense changed-arc recording is 2.0%.

## 2026-07-15 - Reusable changed-node storage

- Replaced the per-solve `std::list` of changed BK nodes with a reusable
  contiguous vector. This removes one allocation and free per changed node;
  all existing consumers only require ordered iteration. The memory estimate
  now includes worst-case retained changed-node and changed-arc indices.
- The 360-solve incremental/full oracle remains exact. Release tests pass in
  32, 64, 128, GMP, legacy replay, and METIS builds, and phase tests pass
  36/36. A zero-loss call profile reduced allocator symbols from roughly 12%
  combined to about 1%; the remaining dominant work is BK itself and exact
  cut-value maintenance.
- Final pinned three-run results against the original implementation are:

| Dataset | Halo | Before (s) | Final (s) | Wall reduction |
| --- | ---: | ---: | ---: | ---: |
| Spiral | h1 | 6.325 | 3.125 | 50.6% |
| Spiral | h2 | 4.484 | 2.506 | 44.1% |
| Spiral | h5 | 10.766 | 6.231 | 42.1% |
| Head | h1 | 12.062 | 5.913 | 51.0% |
| Head | h5 | 11.143 | 6.797 | 39.0% |

- Every before/after run used the same objective, cut count, DD rounds,
  promotions, CPU affinity, and Release binary configuration. h2 remains the
  best Spiral depth; h5 remains over-expanded despite benefiting from the
  common hot-path improvements.

## 2026-07-17 15:13 PDT - Patience-spaced plateau regularization

- Changed all-scale disagreement-plateau regularization from a permanently
  active unit strength to one cumulative unit pulse per complete
  disagreement-patience window. A disagreement improvement resets the window,
  and another pulse cannot occur on the immediately following iteration.
- Kept the strict regularization invariant and existing promotion policy:
  `R >= Q` stops the current scale immediately and promotes the objective by
  10 while retaining the supported warm state.
- Added tracker branch coverage for improvement resets, complete repeated
  windows, and suppression of consecutive pulses. The full Release CTest suite
  passes (2/2).
- A weighted Ghiglia-Pritt Longs direct run with four METIS partitions,
  `Q=2000`, schedule start 125, and disagreement patience 100 certified after
  one promotion, 20,122 DD iterations, and 32 cuts. Full per-iteration progress
  logging raised wall time to 28.1 seconds, so this run validates behavior but
  is not a clean performance result.

## 2026-07-17 15:28 PDT - Configurable scaled-epsilon activation step

- Replaced the fixed `step_size <= 10` check with
  `DualDecompositionOptions::scaled_epsilon_max_step_size`, retaining 10 as
  the default.
- Added direct branch coverage for the default cutoff, a configured cutoff at
  its boundary and immediately above it, the disabled regularization scheme,
  and rejection of nonpositive cutoffs.

## 2026-07-17 16:45 PDT - Configurable scaled-epsilon strength cap

- Added `DualDecompositionOptions::scaled_epsilon_strength_cap`. Zero retains
  the original step-sized strength; a positive value applies
  `min(step_size, cap)` after the configured activation cutoff.
- Kept cumulative budget accounting, strict budget certification, and
  objective-scale promotion unchanged.
- Added branch coverage for active capping, smaller uncapped steps, zero as
  backward-compatible behavior, and negative-value rejection. The native
  partition worker test suite passes.
