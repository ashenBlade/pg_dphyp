# DPHyp join enumeration algorithm for PostgreSQL

Implementation of DP join enumeration algorithm based on hypergraphs for PostgreSQL.

Original paper ["Dynamic Programming Strikes Back"](https://www.researchgate.net/publication/47862092_Dynamic_Programming_Strikes_Back).

## What is it

This is an extension for PostgreSQL that implements dynamic programming algorithm based on hypergraphs to find join ordering (DPhyp).

In PostgreSQL there are 2 builtin algorithms:

- Conventional dynamic programming (called DPsize)
- GEQO, randomized

DPsize can find suitable join ordering, but it's performance dramatically degrades starting from 12 tables.
To prevent infinite awaiting GEQO is used - it is randomized algorithm that is enabled by `geqo` GUC and starts working when tables count hits `geqo_threshold`.

DPhyp algorithm can solve such problem - for some queries it can find *close to optimal* plan in *reasonable amount of time*.

Idea of algorithm is simple: instead of trying to join relations with each other (DPsize) we find links between tables and guide search process using them (DPhyp).
This helps us reduce search space and do not consider obviously impossible JOINs.

> NOTE: found join ordering can be not optimal - DPsize's plan can be more optimal.

## Getting started

Building is the same as for other extensions: clone to `contrib` or use `PGXS` and run `make install`.

Installation to database requires only adding extension to `shared_preload_libraries`:

```text
shared_preload_libraries = 'pg_dphyp'
```

## GUC

There are 2 GUC settings:

1. `pg_dphyp.enabled` (boolean) - controls whether algorithm is enabled. It is `on` by default
2. `pg_dphyp.cj_strategy` (enum) - controls behavior of algorithm related to presence of cross joins:

    - `no` - no actions are taken, if algorithm failed to build final JOIN relation, then DPsize/GEQO is called
    - `pass` - if algorithm failed to build JOIN relation, then it tries to find all built JOIN relations groups (disjoint relation sets) and pass them to DPsize/GEQO (thus we do not throw away done work)
    - `detect` - try to find all disjoint relation sets before algorithm is run and create hyperedges for them. If algorithm still fails to create final JOIN relation, then DPsize/GEQO is called

## Implementation

I assume that you are familiar with DPhyp algorithm, therefore there will be no explanation of it's work.
Instead this section describes key ideas.

### HyperNode representation

PostgreSQL has `Bitmapset` structure that allows us to work with bitmap sets of any size.
But 1) it requires memory allocation each time it is created or modified and 2) in most cases amount of relations in JOIN less than 64.
So for performance reason all operations are performed with `bitmapword`.

To iterate over subsets bit manipulation is used.
Formula `(-x) & x` allows us to iteratively add 1 to existing subset value, so we iterate over incrementing numbers.
i.e. we have neighborhood `010110`, so we will have next subsets:

```text
000010
000100
000110
010000
010010
010100
010110
```

### HyperGraph building

Before DPhyp is run hypergraph must be created and more specifically - hyperedges.

To built hyperedges predicates are used.
List of `RelOptInfo` is passed as input for `join_search_hook` and it has already processed list of `RestrictInfo` that we can use.
Edges are created from predicates found in:

- `RelOptInfo->joinclauses` - generic expressions involving multiple relations
- `PlannerInfo->eclasses` - if there are JOIN clauses involving equality
- `PlannerInfo->join_info_list` - for non-INNER join clauses

But note, that these cases do not cover CROSS JOINs - use `pg_dphyp.cj_strategy` GUC for that (assume that cross joins are not frequent in real workload).

HyperEdges are split into 2 kinds - simple and complex.
Simple hyperedges contains hypernodes with single element, otherwise it is complex.

As an optimization simple edges stored as simple `bitmapword` per relation: in `DPHypContext->simple_edges` array indexed by number of base node and in `HyperNode->simple_edges` - precomputed simple neighborhood for hypernode.

### Cross Join Set

Some join clauses can be complex and contain many different relations which has no explicit join clause between them.
For example:

```sql
t1.x + t2.x = t3.x + t4.x
```

We have hyperedge `{t1, t2} - {t3, t4}`. If there are no join clauses `t1 - t2` and `t3 - t4`, then there is not edge between them in graph.

Such hypernodes in hyperedges called here cross join sets (cjs for short).
All relations in CJS must be 'cross joined' with each other, that is simple hyperedge is created.

This can add complexity and slow down performance, but otherwise we will not be able to create JOIN relation.

### Optimal neighborhood subset iteration

Computing neighborhood of some hypernode is core of the algorithm.
According to [MySQL comments](https://github.com/mysql/mysql-server/blob/ff05628a530696bc6851ba6540ac250c7a059aa7/sql/join_optimizer/subgraph_enumeration.h#L280):

> This function accounts for roughly 20–70% of the total DPhyp running time

So optimization of this function is key to performance.

There are some observations. First, we can compute neighborhood for 1 part of set and then just add some delta for other part.

This allows us to compute neighborhood effectively when some part already known.
In the implementation it is implemented in `get_neighbors_delta` function which accepts some base neighborhood and delta (which must be computed).

Another observation taken from the way we iterate over neighborhood subsets.
Notice, that iteration is just incrementing and when we increment MSB changes much less frequently than LSB, so we can cache that MSB part and calculate only LSB.

The third observation is that to build neighborhood for 1 subset we can use some other subset and just add 1 relation to it, i.e. to build `01010` we can use `00010` or `01000`.

Let's build table of subsets for 4 relations. Arrows will show how we can use previous build subsets to calculate current:

```text
0000 <+  <+  <+
   ^  |   |   |
0001  |   |   |
      |   |   |
0010 -+   |   |
   ^      |   |
0011      |   |
          |   |
0100 <+  -+   |
   ^  |       |
0101  |       |
      |       |
0110 -+       |
   ^          |
0111          |
              |
1000 <+  <+  -+
   ^  |   |
1001  |   |
      |   |
1010 -+   |
   ^      |
1011      |
          |
1100 <+  -+
   ^  |
1101  |
      |
1110 -+
   ^
1111
```

If you look more closely, you will notice this pattern - to calculate current neighborhood we must:

1. Unset first `1` in current subset (i.e. `0101` -> `0100`)
2. Take neighborhood of that subset
3. Process that first `1` with saved 'base' neighborhood

And lastly, note that step 2 is literally number of leading zeros in that number and number of zeros is the number of accesses to it's neighborhood.

But note, that subsets are sparse bitmaps, so instead we use number of current *iteration* - it serves same purposes, but allows comfortable calculation (just represents bitmask of current subset).

So this gives us such caching scheme:

1. Create array of `bitmapword` (sets) of max relation count (in case of PostgreSQL it is `BITS_PER_BITMAPWORD`) - this is our cache.

2. Each time we get new subset remove lowest bit from iteration number, calculate number of zeros and take element from cache at given index - this is base neighborhood
    - If resulting subset is `0` (i.e. when current is `0010`), then take `0` (because when iterating such way old neighborhood is in `excluded` set)
3. Calculate neighborhood for that lowest bit with 'base' neighborhood
4. Calculate number of leading zeros in current iteration number and store resulting neighborhood in cache (at that index)

As you can see this allows us to process single relation at a time and no recalculation required.

### Hyperedge indexing

During calculation of neighborhood we have some moving parts, but other parts are fixed or can grow.

For this reason array of complex edges is indexed - sorted by number of leading zeros in 'right' part of edge.
Example array:

```text
[
    xxx - 0010
    xxx - 1010
    xxx - 0100
    xxx - 1100
    xxx - 1000
]
```

Index is an array where element at `i` points to index in `complex_edges` array which has all relations greater or equal to `i`.
This index is called `start_index`.

For example array we have such index:

```text
[0, 0, 0, 2, 4]
```

We iterate over complex edges in 2 scenarios:

1. Compute neighborhood
2. Decide whether 2 hypernodes are connected

When computing neighborhood `excluded` set is fixed and we known for sure that leading part of it will all be `1`.
So we compute size of this part and start iterating from index where all relations in 'right' part greater or equal to that number.

For example, if `excluded` is `0011`, then index obtained from `start_index` is `2`.

Note that in PostgreSQL we have `bmw_rightmost_one_pos`, but not `bmw_rightmost_zero_pos`, so we use `+-1` trick: add 1 to `excluded` set, calculate rightmost 1 and remove 1 from result value - this gives us number of leading 1.

For second case we don't have `excluded` set, but instead right hypernode is fixed (for which we test for connection).
Idea is the same - 'right' edge will definitely not be subset if it has relations less than minimal from that tested hypernode.

For example, if `right` hypernode is `1000`, then index obtained from `start_index` is `4`.

## Known problems

This algorithm will not find final plan in case of disjoint sets.

They can arise when `CROSS JOIN`/`,` is used without any clauses that connects them.

Also, there is case for outer parameters, because outer parameters are not handled as relations.
For example such subquery has 2 disjoint relations due to outer parameter:

```sql
SELECT * FROM t1
WHERE t1.x IN (SELECT t2.x FROM t2, t3 WHERE t2.x > t1.x AND t3.x <> t1.x);
```
