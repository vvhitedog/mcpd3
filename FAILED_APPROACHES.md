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
