# Failed Approaches

## Capacity precision refactor

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
