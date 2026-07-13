# Capacity Precision Performance

This report measures the runtime and memory cost of the build-time capacity
modes added in July 2026. It also compares current 32-bit behavior with the
32-bit implementation immediately before the precision refactor.

## Conclusion

The current 32-bit path has no material performance regression. Across the
four measured paths, total-wall change versus the pre-refactor checkout ranges
from 2.0% faster to 0.3% slower. A separate 15-pair alternating-order run of
the 64x64 DD case produced a 0.998 median current/baseline ratio.

Wider precision has workload-dependent cost:

- 64-bit capacities add 2.1% to 10.9% total runtime and about 13% peak RSS.
- 128-bit capacities add 17.3% to 317.1% runtime and about 65% to 67% RSS.
- GMP adds 73.1% to 451.8% runtime in the PU tests, 172.9% to 438.1% on the
  Waterloo tests, and about 180% RSS.

The large spread is expected. DD-heavy workloads spend substantial time in
precision-independent scheduling and coordination. Direct BK and graphs whose
maxflow dominates expose arithmetic-width cost more directly.

## Environment

- Date: 2026-07-13 PDT
- CPU: Intel Core i7-9750H, 6 physical cores, 12 hardware threads
- Measured CPUs: physical cores 0 and 1 through `taskset -c 0,1`
- RAM: 15 GiB
- Compiler: GCC 11.4.0
- Build: `Release`, tests disabled in benchmark build directories
- CPU governor observed during the run: `powersave`
- Pre-refactor baseline: `43e0ada`
- Optimized precision implementation: `e6030dd`

The baseline and current sources were checked out into separate detached Git
worktrees, configured into separate build directories, and run with identical
arguments. Mode order was rotated across rounds. Waterloo inputs were read
once before measurement to populate the page cache.

## Inputs And Correctness

The Waterloo inputs use the DIMACS `.max` format and published `.sol` files
from <https://vision.cs.uwaterloo.ca/data/maxflow>.

| Input | Internal nodes | Internal directed arcs | Verified objective |
|---|---:|---:|---:|
| Simulated PU 64x64, seed 1 | 4,096 | 16,002 | 9,383 |
| `LB07-bunny-sml` | 805,800 | 4,782,484 | 961,163 |
| `BL06-gargoyle-sml` | 1,105,920 | 4,276,628 | 29,696,707 |

Every precision mode produced the same labels/objective path and iteration
counts. The PU run used 35 PU iterations, 39 cut attempts, and 8,805 DD
iterations. Bunny and gargoyle used 130 and 131 DD iterations, respectively.

Gargoyle's parser normalization removes a constant terminal imbalance of
17,040,370. The solver reports the normalized objective 12,656,337; adding the
constant gives the published `.sol` value 29,696,707. A one-partition solve
produced the same normalized objective.

Input SHA-256 values:

- `LB07-bunny-sml.max`:
  `2c3f2d18137bbf5fa24dfe2d26feb36ac9bc9129d7f59165cad9c1ed60a5b654`
- `BL06-gargoyle-sml.max`:
  `19448f63126527095734cc4db6c46f493304e3fd0f2ed94d845cc7675e0795ceeff`

## 32-Bit Before And After

Times are medians in milliseconds. Lower is better. Waterloo values use five
runs per build; the PU DD values use five order-rotated runs and were also
checked with 15 alternating current/baseline pairs.

| Workload | Before | Current 32-bit | Change |
|---|---:|---:|---:|
| PU 64x64, mcpd3 native DD | 6,373.250 | 6,306.200 | -1.1% |
| PU 64x64, direct BK | 159.814 | 156.682 | -2.0% |
| Waterloo bunny, total wall | 1,228.896 | 1,232.973 | +0.3% |
| Waterloo gargoyle, total wall | 3,729.835 | 3,727.379 | -0.1% |

The initial precision implementation had a real DIMACS ingestion regression:
it allocated a string and invoked Boost arbitrary-precision parsing for every
capacity, including 32-bit builds. The optimized implementation uses a
single-pass checked native parser and direct checked native conversions.

## Precision Runtime

Times are medians in milliseconds. Parentheses show change from current
32-bit. `Total` includes input, construction, and solve phases. `Solve` excludes
DIMACS reading and graph construction.

| Workload | 32-bit | 64-bit | 128-bit | GMP |
|---|---:|---:|---:|---:|
| PU 64x64, native DD | 6,306.200 | 6,436.120 (+2.1%) | 7,394.980 (+17.3%) | 10,912.900 (+73.1%) |
| PU 64x64, direct BK | 156.682 | 173.707 (+10.9%) | 653.490 (+317.1%) | 864.567 (+451.8%) |
| Bunny, total | 1,232.973 | 1,367.510 (+10.9%) | 2,269.010 (+84.0%) | 6,634.175 (+438.1%) |
| Bunny, solve | 434.788 | 430.830 (-0.9%) | 949.609 (+118.4%) | 3,699.825 (+750.9%) |
| Gargoyle, total | 3,727.379 | 3,908.435 (+4.9%) | 5,305.264 (+42.3%) | 10,173.752 (+172.9%) |
| Gargoyle, solve | 2,926.546 | 2,942.130 (+0.5%) | 3,961.051 (+35.3%) | 7,267.476 (+148.3%) |

The bunny 64-bit solve difference is within run noise; its total cost comes
from wider input and graph construction. GMP construction is expensive because
each capacity is a separately managed arbitrary-precision object.

## Peak RSS

Peak RSS was recorded with `/usr/bin/time`. Values are median MiB; parentheses
show change from current 32-bit.

| Input | 32-bit | 64-bit | 128-bit | GMP |
|---|---:|---:|---:|---:|
| Bunny | 515.3 | 586.8 (+13.9%) | 861.3 (+67.1%) | 1,422.7 (+176.1%) |
| Gargoyle | 493.4 | 556.3 (+12.7%) | 813.1 (+64.8%) | 1,391.5 (+182.0%) |

## Commands

Each precision mode was configured independently:

```bash
cmake -S /path/to/mcpd3 -B build-MODE \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=OFF \
  -DMCPD_CAPACITY_MODE=MODE
cmake --build build-MODE --target mcpd3_native_monolith_benchmark -j2
```

Waterloo command:

```bash
taskset -c 0,1 env MCPD3_PARTITIONER=basic \
  build-MODE/mcpd3_native_monolith_benchmark INPUT.max \
  --directed --partitions 2 --threads 2 \
  --schedule-start 10000 --schedule-levels 5 \
  --max-iterations 10000 --patience 10 --objective-scale 1
```

PU native-DD command:

```bash
taskset -c 0,1 env MCPD3_PARTITIONER=basic \
  build-MODE/phase_unwrapping_benchmark \
  --width 64 --height 64 --seed 1 --anchor 0 \
  --backend mcpd3-n --mcpd-partitions 2 --mcpd3-n-threads 2 \
  --quantum-strategy up-then-down --repeat 1
```

For the direct-BK PU comparison, `--backend mcpd3-n` and its worker arguments
were replaced by `--backend bk`, and each binary performed 15 repeats.
