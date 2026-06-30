# mcpd3

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

## Build And Test

Configure and build:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Run tests:

```bash
ctest --test-dir build --output-on-failure
```

The main test target is `partition_worker_test`, which exercises the extracted
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
