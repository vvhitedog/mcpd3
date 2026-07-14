# Failed Approaches

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

## 2026-07-14 - Undirected PDHG prototype

- Do not include fixed virtual source/sink degrees in the PDHG automatic step
  bound. The primal operator contains only free-vertex columns. Including the
  bunny source's 167,944 terminal edges reduced the step by roughly 155x and
  invalidated the first real-graph timing pass. A regression test now fixes
  the maximum degree to free vertices only.
- Do not enable ergodic primal and dual averages by default on the CPU
  prototype. They add full vertex/edge memory writes every iteration; on the
  16k grid they did not change the 400-iteration certificate and increased
  wall by 36%.
- Do not expect the symmetric unit-scale step to transfer to uniformly large
  capacities. On bunny, reciprocal balance 1000 was far better than balance 1
  or 30000 while preserving `tau*sigma`. Even tuned PDHG remained roughly
  190x slower than BK for exact certification.
