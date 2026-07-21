# Failed Approaches

## Finite promoted schedules

- Do not let a disagreeing unit scale terminate merely because its per-scale
  iteration budget expired. Promote and restart while configured promotion
  headroom remains; only an explicit promotion/resource limit or numeric limit
  may terminate without agreement.
- Do not extend a promoted schedule by a fixed one-level increment. When
  objective scale and initial step differ, promotion can jump multiple
  decades. Recompute the exact number of levels required to reach unit scale.
- Promotion is not a substitute for a usable alpha step. A deliberately poor
  fixed schedule may still reach its explicit promotion limit with primal
  disagreement even when its certified lower bound has reached the objective.


## Persistent quantum state

- Do not treat an arbitrary capacity refresh as a flow-scaling transition.
  Scale local flow only when every internal forward and reverse capacity is
  proven proportional; otherwise preserve/project or reset through the
  existing unscaled path.
- Do not automatically scale DD alpha, last-alpha, or momentum with local
  internal flow. Those variables can include effects from objective terms that
  did not receive the same multiplier.

## Capacity precision refactor

- Do not use the historical 2594-round 128x128 or 3973-round 256x256 quantum
  continuation runs as arithmetic-correct performance baselines. Their
  32-bit `d_flow` accumulation invokes signed overflow, so the wrapped local
  terminal state and derived DD lower bound are not valid C++ or valid solver
  arithmetic even though these instances happened to finish at the known
  primal objective.
- Do not restore 32-bit node balances or Lagrange multipliers to reproduce the
  old fast trajectory. Source capacities and per-arc residuals can remain
  compact, but sums and terminal potentials must use `Objective`.
- Do not raw-copy, `memset`, `realloc`, mmap, or binary-serialize BK nodes,
  arcs, or vectors containing GMP values. GMP objects require construction,
  destruction, and value-aware persistence. Streaming state uses decimal
  length-prefixed integers for nontrivial numeric types.
- Do not use `long` minimum/maximum values as an "unset" objective sentinel.
  GMP has no bounded minimum; explicit presence flags are required.
- Installing `libgmp-dev` system-wide was unavailable in the development
  environment because noninteractive sudo was not authorized. Validation used
  locally extracted distribution packages through `GMP_ROOT` and `BOOST_ROOT`;
  this is a test-environment workaround, not a project dependency layout.
- Do not leave legacy template defaults as `int/long` when the public build has
  a configured precision. The old CSR defaults either narrowed fixed-width
  modes or failed to instantiate in GMP mode.
- Do not raw-mmap capacity-bearing CSR records in GMP mode. The shared array
  API must construct nontrivial values on heap; consequently GMP CSR is exact
  but not an out-of-core capacity store.
- Widening public result typedefs alone is not sufficient. Intermediate
  capacity additions remain signed-overflow hazards, and aggregate flow/cut
  totals must be accumulated in `Objective`; use the checked integer helpers
  at each capacity/objective ownership boundary.
- Do not parse each native DIMACS capacity through a temporary `std::string`
  and `boost::multiprecision::cpp_int`. It made 32-bit DIMACS input roughly 2x
  slower in the first performance pass. Use the checked single-pass native
  token path and reserve arbitrary-precision parsing for wide/GMP values.
- Do not convert native Lagrange updates to capacities through decimal text.
  The conversion is in the DD hot path and caused a repeatable 32-bit runtime
  penalty. Use direct checked numeric conversion.
- Do not infer a small regression from fixed-order benchmark runs. Rotate all
  precision modes across run positions and use alternating baseline/candidate
  pairs when the difference is close to run noise.

## 2026-07-13 - Historical benchmark replay boundary

- Do not enable `MCPD_LEGACY_32BIT_DD_REPLAY` in normal solver builds. It
  exists only to reproduce and audit pre-check benchmark trajectories.
- Do not treat a replayed DD lower bound or final conditioned DOWN objective as
  an arithmetic-valid certificate. Explicit wrapping makes execution defined
  and repeatable, but it does not make the wrapped integer state mathematically
  valid.
- Do not replace the widened production aliases with the replay aliases. The
  replay policy must remain opt-in, 32-bit-only, visibly warned, and covered by
  a normal-build regression test showing that checked overflow remains active.

## 2026-07-14 - Regularization growth and promotion

- Do not represent production scaled-epsilon regularization as a replaceable
  binary anchor. A persistent disagreement can then repeat forever at a fixed
  budget, so neither strict-budget failure nor objective promotion is forced.
- Do not scale BK's opaque residual graph during objective promotion when an
  unscaled regularizer is active. The residual includes flow induced by that
  regularizer, so multiplying it implements `factor * (F + R)` rather than
  `factor * F + R`. Scale explicit primary state and reconstruct the residual
  graph from retained arc flow.
- Do not use the generic `regularized objective - total budget` lower bound as
  the final reported certificate after agreement. It is valid but needlessly
  loose. Agreement plus a strict sub-quantum budget certifies the feasible
  primary objective exactly.

## 2026-07-15 - Halo implementation constraints

- Do not replace `halo_depth=1` with a fully induced one-hop graph. Current
  mcpd3-n uses asymmetric one-owner edge placement; changing that baseline
  would invalidate the required package, trajectory, and performance
  equivalence before deeper halos are evaluated.
- Do not round `capacity / multiplicity`. Compute one checked global LCM
  multiplier and use exact integer factors, or reject the configuration.
- Do not create edge-consensus variables. Duplicate node labels are sufficient
  to make every copied min-cut edge term agree, and edge consensus does not map
  cleanly to the current local primal-dual solver.

## 2026-07-17 - Deferred regularization promotion

- Do not defer objective-scale promotion until an over-budget regularized solve
  reaches agreement. On weighted Ghiglia-Pritt Longs with disagreement
  patience 10, the speculative solve reached the 30,000-iteration limit with
  three disagreements and budget 178,690 against quantum 32,000.
- Choosing a later quantum from the budget observed at speculative agreement,
  including a `Q' >= 2R` headroom rule, does not address cases that stall before
  agreement and adds certificate-state complexity. Preserve immediate strict
  promotion when `R >= Q`.

## 2026-07-17 - Unit alpha cleanup for pinned Shear disagreements

- Do not assume a final no-momentum pass with alpha step 1 resolves the
  generated-coherence Shear pathology. With p4 METIS, objective scale 500,
  step 125, and four promotions, the cleanup run took about 2m14s and stopped
  after 40,000 iterations with the same 84 disagreements.
- A state trace showed no cycle to break: the same 84 global nodes retained the
  same opposing local labels throughout the 30,000-iteration final attempt.
  Momentum advanced alpha by about five units per iteration; exact unit updates
  only traversed the same pinned region more slowly.

## 2026-07-17 - Preserving the pre-promotion alpha-step cap

- Do not preserve `max_step_size=initial_step_size` after multiplying the
  objective and stored alphas during promotion. The unconditional fixed-step
  clamp was an artifact of removing the Polyak policy, not original fixed
  schedule behavior.
- In three-repeat weighted tests, the cap increased Noise from 2,013 to 5,313
  DD rounds and Shear from 5,848 to 32,106 rounds without changing either
  certified objective. Keep the capped branch only for explicit compatibility
  experiments.

## 2026-07-20 - Separate coordinator loop for native streaming

- Do not implement MCPD3-N out-of-core execution by exporting partition
  packages to `PartitionWorkerCoordinator` or by recreating local solvers after
  eviction. That is a separate DD implementation and changed round counts even
  when final objectives matched.
- Do not drop BK residual/search-tree state or primal-dual arc/node flow between
  PU cuts. The accepted resident baseline relies on persistent solver state.
- Native out-of-core mode must execute the existing `DualDecomposition` object
  and differ only in allocation backing. Large persistent arrays, including
  flow and BK residual state, belong in file-backed mappings that remain live.

## 2026-07-20 - Partial worker parity and stale promoted schedules

- Do not treat BK's mmap environment setting as complete worker backing. It
  leaves the enclosing primal-dual topology, capacities, flow, labels, and
  change arrays resident, which is enough to OOM before BK paging can help.
- Do not update only the current objective scale during promotion. Persistent
  PU solves read the configured initial step on the next cut; leaving it stale
  causes native and distributed trajectories to diverge after the first
  promotion even when the current cut happens to finish correctly.
## 2026-07-20 - Product tuning as a low-level default

- Do not change generic `DualDecompositionOptions` objective scale from 1 to
  the product profile's 500. Halo/objective normalization tests and library
  users require neutral defaults. Product tuning belongs in a named entry-point
  policy, shared by MCPD3-N and MCPD4.

## 2026-07-20 - Resident partition-package staging in mmap mode

- Do not build file-backed local graphs and then copy their topology and
  capacities through resident `std::vector` package fields. That creates a
  second partition-sized resident allocation before the worker can map it.
- Do not assume moving a resident package into a file-backed worker applies the
  worker storage policy. The worker must explicitly rehome arrays whose backing
  mode differs, while adopting already mapped arrays without a copy.

## 2026-07-20 - Resident final-label aggregation

- Do not recover a large final cut by asking every partition solve to return a
  resident full-label vector and concatenating those vectors. This creates
  partition-sized worker results plus a duplicate aggregate allocation.
- Recover the already-solved labels through bounded worker copies into the
  coordinator's final mapped destination. The recovery API must not trigger a
  different solve or reconstruct local state.
