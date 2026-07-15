# Halo Dual Decomposition

## Scope

`mcpd3-nh` extends the native dual-decomposition solver with finite graph
halos. `halo_depth=1` is the compatibility mode and must retain the current
`mcpd3-n` package topology, arithmetic, solve trajectory, and performance.

For depths greater than one, every partition owns its original core nodes and
also contains every original node reachable from that core within the selected
number of undirected BFS steps. An original edge appears in every local halo
whose node set contains both endpoints. No independent edge-consensus variable
is introduced; agreement is enforced only between local copies of nodes.

## Exact Objective Accounting

Let `r_v` be the number of partitions containing node `v`, and let `r_e` be
the number containing both endpoints of edge `e`. Define the halo objective
multiplier

```text
Q_h = lcm({r_v} union {r_e}).
```

Every local copy of a node unary receives `Q_h / r_v` times its input unary.
Every local copy of an edge receives `Q_h / r_e` times both directed
capacities. Therefore, when all duplicated labels agree, summing local
objectives gives exactly `Q_h` times the original objective.

The effective objective scale, primary upper bounds, and strict
regularization budget quantum must all include `Q_h`. LCM and capacity
multiplication overflow are hard errors unless an existing explicitly enabled
capacity-saturation policy applies.

## Consensus

Every node present in more than one local halo is a constrained node. The
existing signed pairwise Lagrange constraints are created between its local
copies. Each multiplier appears once with each sign, so all consensus terms
cancel when local objectives are summed. Edges have no separate consensus
terms.

## Compatibility Boundary

The existing `mcpd3-n` one-hop representation assigns each edge to one owner
partition and creates only the remote endpoint clone required by that edge.
That is retained verbatim for `halo_depth=1`, even though it is more compact
than constructing a fully induced one-hop halo. Full BFS halo construction and
multiplicity weighting begin at depth two.

## Verification Gates

1. `halo_depth=1` package exports and native solve snapshots are identical to
   the pre-halo implementation.
2. Tiny finite and infinite halo layouts have exact node/edge memberships and
   checked objective multipliers.
3. Summed local objectives equal `Q_h` times direct min-cut objectives on
   exhaustive and randomized graphs.
4. Capacity replacement, objective promotion, flow-heat reconstruction,
   package-only export, and worker execution retain their existing contracts.
5. The 64x64, 128x128, and 256x256 PU benchmarks show no material h1 runtime
   regression before any h2+ performance conclusion is accepted.
