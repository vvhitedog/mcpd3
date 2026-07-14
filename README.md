# mcpd3

## Reference-Guided Local Mincuts

`DualDecompositionOptions::reference_cut_labels` optionally supplies one
binary global reference label per node. Each partition maps those labels to its
local copies and uses them only to select an exact local mincut.

`ReferenceCutSelection::EXACT_REFERENCE_IF_OPTIMAL` accepts the reference when
it satisfies the solved residual graph and otherwise keeps BK's exact cut.
`ReferenceCutSelection::CLOSEST_EXACT` instead solves a residual minimum-closure
problem to find the exact local mincut with minimum Hamming distance to the
reference. `reference_cut_check_interval` controls how often the selection is
applied and must be positive.

Because both modes preserve every local Lagrangian minimum value, the summed DD
lower bound remains valid. Primal agreement therefore retains the normal exact
binary-mincut certificate.

mcpd3 is a C++ min-cut/max-flow solver library and benchmark executable set.
The name is pronounced "mcpd cubed": minimum cut, primal-dual, dual
decomposition.

The repository can be used in two ways:

- As a standalone min-cut solver for DIMACS max-flow/min-cut graphs.
- As the solver core for a distributed runtime such as
  [mcpd4](https://github.com/vvhitedog/mcpd4), which handles TCP
  coordinator/worker process orchestration while mcpd3 owns the optimization
  logic.

## What It Does

mcpd3 solves s-t min-cut problems. It includes:

- a direct `PrimalDualMinCutSolver`;
- DIMACS graph readers, including streaming readers for larger inputs;
- graph partitioning helpers;
- a dual-decomposition optimizer that splits a graph into partition
  subproblems and coordinates boundary-node agreement;
- a partition-worker API used by distributed runtimes;
- exact scaled-epsilon regularization diagnostics and objective-scale
  promotion for regularized agreement recovery.

The implementation uses the Boykov-Kolmogorov maxflow code under `maxflow/`
for local subproblems and wraps it with mcpd3 graph, primal-dual, and
dual-decomposition abstractions.

## Repository Layout

- `primaldual/`: direct primal-dual min-cut solver API.
- `graph/`: graph containers, DIMACS readers, and partitioning.
- `decomp/`: dual decomposition, lower-bound accounting, and partition-worker
  interfaces.
- `maxflow/`: local maxflow implementation used by the solvers.
- `example/`: standalone executables.
- `tests/`: CTest-based coverage for partition workers, coordinator behavior,
  regularization, and objective-scale promotion.

## Experimental Undirected PDHG

The `exp/pdhg-undirected-mincut` branch contains a local CPU prototype of PDHG
for the exact convex relaxation of undirected s-t min-cut. It is optional and
does not replace the BK-backed `PrimalDualMinCutSolver`.

`MinCutGraph` inputs must have equal forward and reverse capacity for every
internal edge. Signed terminal capacities are expanded to fixed virtual source
and sink edges. Directed/asymmetric inputs are rejected rather than
symmetrized.

Build and compare both solvers:

```bash
cmake -S . -B build-pdhg -DCMAKE_BUILD_TYPE=Release \
  -DMCPD_BOOST_INCLUDE_DIR=/home/matt/.local/boost-dev/usr/include
cmake --build build-pdhg -j

./build-pdhg/mcpd3_pdhg_benchmark \
  --dimacs graph.max \
  --solver both \
  --pdhg-check-interval 100 \
  --pdhg-max-iterations 100000
```

`--solver existing`, `--solver pdhg`, and `--solver both` select the
BK-backed solver, PDHG, or comparison mode. Important experimental controls
include:

```text
--pdhg-tau X --pdhg-sigma X --pdhg-theta X
--pdhg-step-size-scale X --pdhg-step-balance X
--pdhg-use-ergodic-primal 0|1 --pdhg-use-ergodic-dual 0|1
--pdhg-max-iterations N --pdhg-check-interval N
--pdhg-time-limit-seconds X
--pdhg-stagnation-checks N --pdhg-stagnation-tolerance X
--capacity-quantum X --pdhg-eps-abs X --pdhg-eps-rel X
--pdhg-lower-bound-safety-factor X
--print-history 0|1
```

The default automatic steps are
`tau=sigma=0.99/sqrt(2*maximum_free_vertex_degree)`. `step_balance=b` changes
them to `tau/b` and `sigma*b`, preserving their product. This is useful when
capacities are uniformly much larger than one. Explicit `tau` or `sigma`
overrides the corresponding automatic value.

Every reported cut uses exact configured integer capacities. The dual lower
bound uses directed rounding plus a conservative numerical margin. Exactness
is reported only when the certified gap is strictly below the configured
capacity quantum. Iteration, time, and stagnation limits remain explicitly
uncertified statuses.

The current CPU results are experimental. A 16,384-node grid certifies in
0.297 seconds total versus 0.0075 seconds for BK. An 805,800-node bunny
instance with `--pdhg-step-balance 1000` finds the exact cut after 27.68
seconds of solving and certifies after 32.68 seconds, versus 0.168 seconds of
BK solve time. PDHG is therefore not currently competitive with BK as a local
exact solver.

## Build And Test

Configure and build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Capacity precision is selected at configure time with
`MCPD_CAPACITY_MODE`. The default remains `32` for compatibility:

```bash
# 32-bit capacities, 64-bit accumulated objectives (default)
cmake -S . -B build-32 -DMCPD_CAPACITY_MODE=32

# 64-bit capacities, 128-bit accumulated objectives
cmake -S . -B build-64 -DMCPD_CAPACITY_MODE=64

# 128-bit capacities, 256-bit accumulated objectives
cmake -S . -B build-128 -DMCPD_CAPACITY_MODE=128

# Arbitrary-precision capacities and objectives using GMP
cmake -S . -B build-gmp -DMCPD_CAPACITY_MODE=gmp
```

All modes require Boost.Multiprecision headers. GMP mode additionally requires
the GMP C and C++ development packages, commonly installed as `libgmp-dev` on
Debian/Ubuntu. `BOOST_ROOT` and `GMP_ROOT` can point CMake at non-system
installations.

`Capacity` is the compact source-data and arc-residual type selected above.
Accumulated objectives, per-node flow balances, terminal residuals, and DD
Lagrange multipliers use the next wider `Objective` type. This keeps 32-bit
arc storage compact while preventing valid sums and optimization state from
being narrowed back into 32 bits.

GMP capacities are nontrivial C++ objects, so BK node/arc arrays and
capacity-bearing CSR arrays use constructed heap storage in GMP mode.
File-backed and anonymous mmap storage remain available for the fixed-width
modes. The GMP CSR path is exact but is not an out-of-core capacity store.

See [CAPACITY_PERFORMANCE.md](CAPACITY_PERFORMANCE.md) for measured runtime and
memory costs of each mode, including a pre/post-refactor 32-bit comparison.

Run tests:

```bash
ctest --test-dir build --output-on-failure
```

`capacity_precision_test` drives the configured extreme value through BK,
DIMACS parsing, the mcpd3 solver, in-process workers, and streaming worker
storage. `partition_worker_test` exercises the extracted
partition-worker API, in-process coordinator loop, batched worker solves,
scaled-epsilon regularization, objective-scale promotion, randomized initial
alpha behavior, and lower-bound certificate accounting.

Optional build flags:

```bash
cmake -S . -B build -DMETIS_ENABLED=ON
cmake -S . -B build -DTRAP_SIGNED_INTEGER_OVERFLOW=ON
cmake -S . -B build -DGPERF_PROFILER_BUILD=ON
```

`METIS_ENABLED=ON` enables METIS partitioning if METIS and GKlib are
available. Without METIS, mcpd3 falls back to the built-in basic/local
partitioning paths.

## Standalone Examples

Run the built-in tiny graph example:

```bash
./build/simple_example
```

Run direct DIMACS solving:

```bash
./build/dimacs_example /path/to/graph.max
```

Run the dual-decomposition DIMACS driver:

```bash
./build/dimacs_dual_decomp_example /path/to/graph.max \
  --partitions 10 \
  --max-iterations 10000 \
  --threads 4 \
  --capacity-multiplier 10000 \
  --disable-primal-upper-bound
```

For directed DIMACS inputs, use the directed streaming reader:

```bash
./build/dimacs_dual_decomp_example /path/to/graph.max \
  --stream-directed-input \
  --partitions 10 \
  --capacity-multiplier 10000
```

Useful dual-decomposition options:

```text
--partitions N
--patience N
--max-iterations N
--threads N
--regularization scaled-epsilon|none
--regularization-budget-limit N
--disable-scale-promotion
--max-scale-promotions N
--random-initial-alpha-radius N
--random-initial-alpha-seed N
--capacity-multiplier N
--stream-directed-input
--stream-symmetric-input
--disable-primal-upper-bound
--quiet
```

## Native Dual-Decomposition Benchmark

`mcpd3_native_monolith_benchmark` is the clean native-local comparator for
dual decomposition. It constructs `mcpd3::DualDecomposition` directly and, by
default, disables partition-package export so local benchmarking does not copy
subproblems into the distributed worker representation.

Build it with the normal mcpd3 build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Run a directed DIMACS benchmark:

```bash
MCPD3_PARTITIONER=basic ./build/mcpd3_native_monolith_benchmark \
  /data/adhead.n6c10.max \
  --directed \
  --partitions 10 \
  --objective-scale 1000 \
  --schedule-start 10000 \
  --schedule-levels 5 \
  --max-iterations 10000 \
  --exhaust-regularized-scale-iterations
```

The benchmark prints unbuffered key/value fields for graph size, memory
snapshots, objective values, agreement status, objective-scale promotions, and
separate read/scale/construct/solve timings. Use `--emit-partition-packages`
only for diagnostics that intentionally compare against the worker export path.
When comparing to mcpd4 defaults, pass
`--exhaust-regularized-scale-iterations` so the low-scale schedule matches the
worker-coordinator path.
`--saturate-capacity-overflow` and `--truncate-capacity-overflow` are opt-in
compatibility modes that clamp scaled capacities to the 32-bit range, including
during later objective-scale promotions.

## Programmatic Use

For direct in-process solving, include the repository root and link the maxflow
sources:

```cpp
#include <graph/dimacs.h>
#include <primaldual/mcpd3.h>

int main() {
  auto graph = mcpd3::read_dimacs("graph.max");
  mcpd3::PrimalDualMinCutSolver solver(std::move(graph));
  solver.solve();
  return solver.getMinCutValue() < 0;
}
```

For dual decomposition:

```cpp
#include <decomp/dualdecomp.h>
#include <graph/dimacs.h>

int main() {
  auto graph = mcpd3::read_dimacs("graph.max");
  mcpd3::DualDecompositionOptions options;
  options.thread_count = 4;
  options.track_primal_upper_bound = false;
  options.emit_partition_packages = false;
  options.objective_scale = 10000;

  mcpd3::DualDecomposition solver(/*npartition=*/10, std::move(graph),
                                  options);
  solver.solve<true>([](const std::vector<bool> &, double,
                        const std::list<int> &) { return false; });
  return solver.getLastDisagreementCount() == 0 ? 0 : 1;
}
```

When embedding mcpd3 directly, compile/link `maxflow/graph.cpp` and
`maxflow/maxflow.cpp` with your target and add the repo root to the include
path.

## Partition-Worker API

The productized branch exposes a network-free worker API in
`decomp/partition_worker.h` and `decomp/partition_coordinator.h`.

Important types:

- `mcpd3::PartitionPackage`: serialized local subproblem data plus boundary
  constraint endpoints.
- `mcpd3::PartitionWorker`: abstract worker interface.
- `mcpd3::InProcessPartitionWorker`: local implementation of the worker
  interface.
- `mcpd3::PartitionWorkerCoordinator`: coordinator loop that owns alpha state,
  dispatches partition solve rounds, gathers labels/lower bounds, updates
  multipliers, handles regularization, and performs objective-scale promotion.
- `mcpd3::PartitionWorkerCoordinatorOptions`: solver schedule and diagnostic
  options for the worker-coordinator path.

`DualDecomposition::getPartitionPackages()` exports the partition packages
needed to run the same subproblems through the worker API.

## Relationship To mcpd4

mcpd3 deliberately does not own TCP, process management, deployment, or
wire-format concerns. It provides the solver core and network-free partition
worker contracts.

mcpd4 builds on that boundary:

- mcpd3 reads/partitions graphs and defines partition solve requests/results.
- mcpd3 provides the in-process coordinator and worker semantics used as the
  correctness reference.
- mcpd4 serializes mcpd3 partition packages and solve messages over TCP.
- mcpd4 owns coordinator/worker binaries, worker handshakes, batching, process
  integration tests, and runbooks.

This split keeps mcpd3 usable as an independent solver while allowing mcpd4 to
productize distributed execution without adding networking dependencies to the
solver repo.

## Runtime Notes

- `MCPD3_PROGRESS=1` enables partitioning progress reports.
- `MCPD3_PARTITIONER=basic|contiguous|local` selects the built-in partitioner
  mode when METIS is not used.
- `MCPD3_LOCAL_PARTITION_PASSES`,
  `MCPD3_LOCAL_PARTITION_LAMBDA`, and
  `MCPD3_LOCAL_PARTITION_BALANCE_SLACK` tune the local partitioner.
- `--capacity-multiplier` is also the objective scale used by exact
  scaled-epsilon regularization. Larger values give more regularization
  resolution but increase 32-bit capacity overflow risk.
- `--disable-primal-upper-bound` is useful for lower-bound/dual-decomposition
  benchmarking when primal decoding is not needed.

## License

mcpd3 is distributed under the GNU General Public License. See
[LICENSE](LICENSE) and the bundled maxflow license files under `maxflow/`.
