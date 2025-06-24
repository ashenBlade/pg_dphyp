#include "postgres.h"

#include <limits.h>

#include "fmgr.h"
#include "funcapi.h"
#include "miscadmin.h"
#include "nodes/bitmapset.h"
#include "optimizer/geqo.h"
#include "optimizer/pathnode.h"
#include "optimizer/paths.h"
#include "utils/builtins.h"
#include "utils/guc.h"
#include "utils/hsearch.h"

#include "simplebms.h"

#if PG_MAJORVERSION_NUM < 13
#define generate_useful_gather_paths(x,y,z) generate_gather_paths(z,y,z)
#endif

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(pg_dphyp_get_statistics);
PG_FUNCTION_INFO_V1(pg_dphyp_reset_statistics);

void _PG_init(void);
void _PG_fini(void);

/*
 * Represents a hypernode in query graph
 */
typedef struct HyperNode
{
	/*
	 * Bitmap of relations this Hypernode represents.
	 * Must be first field in structure - used as key.
	 */
	bitmapword set;

	/*
	 * Cached representative node of this Hypernode as used in original paper.
	 */
	int representative;

	/*
	 * Created 'RelOptInfo' for this HyperNode.
	 * During DPhyp algorithm this is not NULL only for base hypernodes.
	 * At the end of algorithm we build RelOptInfo for all
	 */
	RelOptInfo *rel;

	/*
	 * List of Hypernode pairs, that can contribute for creating this HyperNode.
	 * Used as an indicator that this HyperNode has a plan and can be created,
	 * even if actually 'make_join_rel' will not be able to create RelOptInfo
	 * from it.
	 */
	List *candidates;

	/*
	 * Cached bitmap of nodes that are connected with this HyperNode
	 * with simple edges.
	 * Just bit OR of 'simple_edges' of all 'nodes'.
	 */
	bitmapword simple_edges;
} HyperNode;

/*
 * Pair of 'left' and 'right' bitmapwords representing hypernodes.
 * lhs and rhs must not intersect.
 * Each edge stored twice - with lhs/rhs swapped, so no need to swap
 * every time and check again.
 */
typedef struct HyperEdge
{
	/* Left side of edge */
	bitmapword left;
	/* Right side of edge */
	bitmapword right;
} HyperEdge;

/*
 * Array of Hyperedges stored as plain C-array instead of List *.
 * Entries stored sorted by tuple (left, right) edges.
 */
typedef struct EdgeArray
{
	/* Capacity of 'edges' array */
	int capacity;
	/* Actual size of 'edges' array (payload) */
	int size;
	/* Array of hyperedges */
	HyperEdge *edges;
	/* Size of 'start_idx' array */
	int start_idx_size;
	/* Index storing positions from which to start iterating */
	int8 *start_idx;
} EdgeArray;

/*
 * Context object, that passed along with any function invocation.
 */
typedef struct DPHypContext
{
	/*
	 * Original planner info
	 */
	PlannerInfo *root;

	/*
	 * List of initial passed 'RelOptInfo' objects.
	 */
	List *initial_rels;

	/*
	 * List of Hypernodes created for initial relations.
	 */
	List *base_hypernodes;

	/*
	 * Size of 'simple_edges' and 'complex_edges' arrays
	 */
	int edges_size;

	/*
	 * Map of base hypernode to bitmaps containing all nodes
	 * that node has simple edge with.
	 */
	bitmapword *simple_edges;

	/*
	 * Array of hyperedges appearing in query graph. Index is id
	 * of hypernode and each array contains all hyperedges that
	 * this node appears in.
	 */
	EdgeArray *complex_edges;

	/*
	 * Dynamic programming table that maps: bitmapword -> HyperNode.
	 */
	HTAB *dptable;
} DPHypContext;

/*
 * Structure used as state for enumerating subsets of given bitmap
 */
typedef struct SubsetIteratorState
{
	/* Common state for subset iteration */
	/*
	 * Current subset value
	 */
	bitmapword subset;
	/*
	 * Current subset to return. 0 means no more subsets.
	 */
	bitmapword state;
	/*
	 * Initial bitmap that used as mask to iterate.
	 */
	bitmapword init;

	/* State to compute neighborhood using deltas */
	/*
	 * Current iteration number. Used to decide which actions to take.
	 */
	bitmapword iteration;

	/*
	 * Cache of neighborhoods for different subsets.
	 * Indexed by number of leading zeros in subset.
	 */
	bitmapword cached_neighborhood[BITS_PER_BITMAPWORD];
} SubsetIteratorState;

#define IS_ODD(number) (((number) & 1) == 1)

/* Collected statistics during backend work */
typedef struct Statistics
{
	/* Total amount of DPhyp was executed */
	uint64 total_runs;
	/* Amount of attempts to perform join search */
	uint64 failed_runs;
} Statistics;

typedef enum CrossJoinStrategy
{
	/* No actions are performed - transfer control to DPsize/GEQO */
	CJ_STRATEGY_NO,
	/* 
	 * Perform check during edges initialization and create hyperedges for
	 * disjoint hypernodes (relations).
	 */
	CJ_STRATEGY_DETECT,
	/* Collect all relations we were able to create and pass them to DPsize/GEQO */
	CJ_STRATEGY_PASS,
} CrossJoinStrategy;

static const struct config_enum_entry cross_join_strategy_options[] =
{
	{ "no", CJ_STRATEGY_NO, false },
	{ "detect", CJ_STRATEGY_DETECT, false },
	{ "pass", CJ_STRATEGY_PASS, false },
	{ NULL, 0, false },
};

/* GUC */
/* Extension is enabled and should run DPhyp */
static bool dphyp_enabled = true;
/*
 * In case of CROSS JOINs we can get disjoint subgraphs for tree, so let user
 * decide how to handle them.
 */
static int dphyp_cj_strategy = CJ_STRATEGY_PASS;

static join_search_hook_type prev_join_search_hook = NULL;

static Statistics statistics = {0};

/* Initialization */
static bitmapword map_to_internal_bms(List *initial_rels, Bitmapset *original);
static void initialize_edges(PlannerInfo *root, List *initial_rels,
							 DPHypContext *context);
static void distribute_cjs(DPHypContext *context, bitmapword cjs);
static void distribute_hyperedge(DPHypContext *context, HyperEdge edge);
static void initialize_hypernodes(DPHypContext *context);
static void hyperedge_array_add(EdgeArray *array, HyperEdge edge);
static void compute_start_index(DPHypContext *context);
static HyperEdge hyperedge_swap(HyperEdge edge);

/* Runtime */
static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes);
static void subset_iterator_init(SubsetIteratorState *state, bitmapword bmw);
static bool subset_iterator_next(SubsetIteratorState *state);
static bitmapword get_neighbors_iter(DPHypContext *context, bitmapword subgroup, bitmapword excluded,
									  SubsetIteratorState *iter_state);
static bitmapword get_neighbors(DPHypContext *context, HyperNode *node, bitmapword excluded);
static int get_start_index(EdgeArray *edges, bitmapword bmw);
static bool hypernode_has_direct_edge_with(DPHypContext *context, HyperNode *node, int id);
static bool hypernode_has_edge_with(DPHypContext *context, HyperNode *node, bitmapword bms);
static void emit_csg_cmp(DPHypContext *context, HyperNode *subgroup,
						 HyperNode *complement);
static void enumerate_cmp_recursive(DPHypContext *context, HyperNode *node,
									HyperNode *complement, bitmapword excluded,
									bitmapword neighborhood);
static void emit_csg(DPHypContext *context, HyperNode *node, bitmapword neighborhood);
static void enumerate_csg_recursive(DPHypContext *context, HyperNode *node,
									bitmapword excluded, bitmapword neighborhood);
static void solve(DPHypContext *context);
static RelOptInfo *hypernode_get_rel(DPHypContext *context, HyperNode *node);

/* Extension */
static RelOptInfo *dphyp(DPHypContext *context, PlannerInfo *root, List *initial_rels);
static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed,
									 List *initial_rels);

Datum
pg_dphyp_get_statistics(PG_FUNCTION_ARGS)
{
	TupleDesc tupdesc;
	Datum values[2] = {0};
	bool nulls[2] = {0};

	if (get_call_result_type(fcinfo, NULL, &tupdesc) != TYPEFUNC_COMPOSITE)
		elog(ERROR, "return type must be a row type");

	values[0] = Int64GetDatum(statistics.total_runs);
	values[1] = Int64GetDatum(statistics.failed_runs);

	PG_RETURN_DATUM(HeapTupleGetDatum(heap_form_tuple(tupdesc, values, nulls)));
}

Datum
pg_dphyp_reset_statistics(PG_FUNCTION_ARGS)
{
	statistics.failed_runs = 0;
	statistics.total_runs = 0;
	PG_RETURN_VOID();
}

static inline bool
hyperedge_is_simple(HyperEdge edge)
{
	return bmw_is_singleton(edge.left) && bmw_is_singleton(edge.right);
}

static inline bool
hyperedge_is_valid(HyperEdge edge)
{
	/*
	 * Vertexes must be not empty and they must not intersect
	 */
	return !(   bmw_is_empty(edge.left)
			 || bmw_is_empty(edge.right)
			 || bmw_overlap(edge.left, edge.right));
}

static inline int
hyperedge_cmp(HyperEdge a, HyperEdge b)
{
	/* Simple integer tuple (lowest(right), left, right) comparison */
	bitmapword t;

	/* Use lowest_bit instead of bmw_first - same result, but faster */
	t = bmw_lowest_bit(a.right) - bmw_lowest_bit(b.right);
	if (t != 0)
		return t;

	t = a.left - b.left;
	if (t != 0)
		return t;

	t = a.right - b.right;
	return t;
}

/*
 * Check that we calculated any query plan for this hypernode
 */
static inline bool
hypernode_has_rel(HyperNode *node)
{
	return node->rel != NULL || node->candidates != NIL;
}

/*
 * Get neighborhood of given hypernode excluding 'excluded' set.
 * Neighborhood calculated from scratch.
 */
static bitmapword
get_neighbors(DPHypContext *context, HyperNode *node, bitmapword excluded)
{
	bitmapword neighbors;
	int idx;

	excluded |= node->set;
	neighbors = node->simple_edges;

	idx = -1;
	while ((idx = bmw_next_member(node->set, idx)) >= 0)
	{
		EdgeArray *complex_edges = &context->complex_edges[idx];
		if (0 < complex_edges->size)
		{
			int i = get_start_index(complex_edges, neighbors | excluded);
			for (; i < complex_edges->size; i++)
			{
				HyperEdge edge = complex_edges->edges[i];
				if ( bmw_is_subset(edge.left, node->set) &&
					!bmw_overlap(edge.right, neighbors | excluded))
				{
					neighbors |= bmw_lowest_bit(edge.right);
				}
			}
		}
	}

	neighbors = bmw_difference(neighbors, excluded);

	return neighbors;
}

static inline bitmapword
get_parent_neighborhood(DPHypContext *context, SubsetIteratorState *iter_state)
{
	int zero_count;
	bitmapword last_bit_removed;
	
	Assert(iter_state->iteration != 0);
	last_bit_removed = bmw_difference(iter_state->iteration, bmw_lowest_bit(iter_state->iteration));
	if (unlikely(bmw_is_empty(last_bit_removed)))
	{
		/* Don't have any parent */
		return 0;
	}

	zero_count = bmw_rightmost_one_pos(last_bit_removed);
	return iter_state->cached_neighborhood[zero_count];
}

/*
 * Get bitmap of neighbors for node excluding all specified.
 * Corresponds to 'N(S, X)' function in paper.
 */
static bitmapword
get_neighbors_iter(DPHypContext *context, bitmapword subgroup,
					bitmapword excluded, SubsetIteratorState *iter_state)
{
	int idx;
	bitmapword neighbors;
	EdgeArray *complex_edges;

	excluded |= subgroup;

	iter_state->iteration++;
	Assert(!bmw_is_empty(iter_state->subset));
	idx = bmw_rightmost_one_pos(iter_state->subset);

	/* Starting point is parent neighborhood */
	neighbors = get_parent_neighborhood(context, iter_state);

	/* Add simple neighborhood */
	neighbors |= bmw_difference(context->simple_edges[idx], excluded);

	/* And neighbors from complex edges */
	complex_edges = &context->complex_edges[idx];

	if (0 < complex_edges->size)
	{
		int i = get_start_index(complex_edges, neighbors | excluded);
		for (; i < complex_edges->size; i++)
		{
			HyperEdge edge = complex_edges->edges[i];
			if ( bmw_is_subset(edge.left, subgroup) &&
				!bmw_overlap(edge.right, neighbors | excluded))
			{
				neighbors |= bmw_lowest_bit(edge.right);
			}
		}

	}
	
	neighbors = bmw_difference(neighbors, excluded);

	/* 
	 * Save current neighborhood in table, but for performance
	 * skip odd-numbered iterations, because they are not used anymore.
	 */
	if (!IS_ODD(iter_state->iteration))
	{
		int zero_count;

		zero_count = bmw_rightmost_one_pos(iter_state->iteration);
		iter_state->cached_neighborhood[zero_count] = neighbors;
	}


	return neighbors;
}

/*
 * Check that 'node' has direct edge with node 'id'.
 * This is not the same as 'has_edge_with' because we must check
 * that it has simple edge
 */
static bool
hypernode_has_direct_edge_with(DPHypContext *context, HyperNode *node, int id)
{
	EdgeArray *edges;
	bitmapword right_bmw;

	/* If we have direct simple edge, then we are done */
	if (bmw_is_member(node->simple_edges, id))
		return true;

	/* Otherwise, we may have complex edge with single 'id' node at right side */
	edges = &context->complex_edges[id];

	right_bmw = bmw_make_singleton(id);

	/* For clarification see comment in 'hypernode_has_edge_with' */
	if (0 < edges->size)
	{
		int i = get_start_index(edges, node->set - 1);
	
		for (; i < edges->size; i++)
		{
			HyperEdge edge = edges->edges[i];
			if (edge.left != right_bmw)
				continue;
	
			if (bmw_is_subset(edge.right, node->set))
				return true;
		}
	}

	return false;
}

/*
 * Check that 'node' has any edge that can be used as connection to 'bmw'.
 * This is used to check that subgroup and complement can be connected
 * to further call 'emit_csg_cmp' and create join rel for them.
 */
static bool
hypernode_has_edge_with(DPHypContext *context, HyperNode *node, bitmapword bmw)
{
	int idx;

	Assert(!bmw_overlap(node->set, bmw));

	/* Check that we have simple edges that connect to 'bmw' */
	if (bmw_overlap(node->simple_edges, bmw))
		return true;

	/* Now check any complex edge has connection to 'bmw' */
	idx = -1;
	while ((idx = bmw_next_member(node->set, idx)) >= 0)
	{
		EdgeArray *edges = &context->complex_edges[idx];

		/* 
		 * We should find edge 'right' side of which is subset of 'bmw'.
		 * You can see that we can skip all 'right's if it has any relations
		 * less than any relation from 'bmw'. Simple case:
		 * 
		 * right: 000110
		 *   bmw: 001100
		 * 
		 * 'right' is not subset of 'bmw' because of second bit set, so
		 * using 'start_index' on 'right' we can quickly find start index
		 * for traverse.
		 * 'get_start_index' works with 'excluded', but not 'bmw' (as in this
		 * case), but as you can see all leading zeros in 'bmw' can be used
		 * as excluded bitmask - we just decrement 1 and all leading 0 become 1:
		 * 001100 -> 001011
		 */
		if (0 < edges->size)
		{
			int i = get_start_index(edges, bmw - 1);
			for (; i < edges->size; i++)
			{
				HyperEdge edge = edges->edges[i];
				if (bmw_is_subset(edge.left, node->set) &&
					bmw_is_subset(edge.right, bmw))
					return true;
			}
		}
	}

	return false;
}

static void
subset_iterator_init(SubsetIteratorState *state, bitmapword base_neighborhood)
{
	state->init = base_neighborhood;
	state->state = (-base_neighborhood) & base_neighborhood;
	state->subset = 0;
	state->iteration = 0;
}

static bool
subset_iterator_next(SubsetIteratorState *state)
{
	if (state->state == 0)
		return false;

	state->subset = state->state;
	state->state = (state->state - state->init) & state->init;
	return true;
}

/* Store 'subgraph'/'complement' pair to further use them to search query plan  */
static void
emit_csg_cmp(DPHypContext *context, HyperNode *subgraph, HyperNode *complement)
{
	HyperNode *hypernode;

	/*
	 * Now we do not create 'RelOptInfo' for this join, but instead
	 * save pair of hypernodes that can be joined together.
	 *
	 * PostgreSQL's planner designed highly cohesion with DPsize algorithm,
	 * so during processing 1 level of join we just call 'make_join_rel'
	 * with nodes of lower level and add more available paths and at the
	 * end we call 'set_cheapest' to find best paths among discovered.
	 * It would be easier to code to just call 'make_join_rel' here and
	 * 'set_cheapest' at the end, but we can not do this, because 'make_join_rel'
	 * expects that 'set_cheapest' was already called with rel at lower level.
	 * So adding 'make_join_rel' + 'set_cheapest' (and some other functions)
	 * here will add overhead by calling them multiple times for same rel.
	 */
	hypernode = get_hypernode(context, subgraph->set | complement->set);
	if (hypernode->candidates != NIL)
	{
		hypernode->candidates = lappend(hypernode->candidates, subgraph);
		hypernode->candidates = lappend(hypernode->candidates, complement);
	}
	else
	{
		hypernode->candidates = list_make2(subgraph, complement);
	}
}

/* 
 * For given 'complement' of 'subgraph' try to enlarge 'complement' using
 * it's neighborhood.
 */
static void
enumerate_cmp_recursive(DPHypContext *context, HyperNode *subgraph, HyperNode *complement,
						bitmapword excluded, bitmapword neighborhood)
{
	SubsetIteratorState subset_iter;

	Assert(!bmw_is_empty(neighborhood));

	subset_iterator_init(&subset_iter, neighborhood);
	while (subset_iterator_next(&subset_iter))
	{
		HyperNode *expanded_complement;

		expanded_complement = get_hypernode(context, complement->set | subset_iter.subset);

		if (hypernode_has_rel(expanded_complement) &&
			hypernode_has_edge_with(context, subgraph, expanded_complement->set))
			emit_csg_cmp(context, subgraph, expanded_complement);
	}

	excluded |= neighborhood;

	subset_iterator_init(&subset_iter, neighborhood);
	while (subset_iterator_next(&subset_iter))
	{
		HyperNode *expanded_complement;
		bitmapword current_neighborhood;
		current_neighborhood = get_neighbors_iter(context, complement->set | subset_iter.subset,
											 	  excluded, &subset_iter);
		if (bmw_is_empty(current_neighborhood))
			continue;

		expanded_complement = get_hypernode(context, complement->set | subset_iter.subset);
		enumerate_cmp_recursive(context, subgraph, expanded_complement,
								excluded, current_neighborhood);
	}
}

/* Find complement for specified 'subgraph' */
static void
emit_csg(DPHypContext *context, HyperNode *subgraph, bitmapword neighborhood)
{
	bitmapword excluded;
	int i;

	Assert(!bmw_is_empty(neighborhood));

	excluded = subgraph->set | bmw_all_bit_set(subgraph->representative);

	i = -1;
	while ((i = bmw_prev_member(neighborhood, i)) >= 0)
	{
		HyperNode *complement;
		bitmapword excluded_ext;
		bitmapword complement_neighborhood;

		complement = (HyperNode *) list_nth(context->base_hypernodes, i);

		/*
		 * Here in original paper we create S = {v} and then check that
		 * edge rhs is subset of S.  But as you can see subset of single element
		 * set is that set itself, so we can make optimized searching
		 * for such edge.
		 */
		if (hypernode_has_direct_edge_with(context, subgraph, i))
			emit_csg_cmp(context, subgraph, complement);

		/*
		 * We are iterating backwards on neighbors, so we have to exclude
		 * all nodes lower current, otherwise, we will get duplicates
		 * and execution time will skyrocket.
		 */
		excluded_ext = excluded | (neighborhood & bmw_all_bit_set(i));
		complement_neighborhood = get_neighbors(context, complement, excluded_ext);
		if (!bmw_is_empty(complement_neighborhood))
			enumerate_cmp_recursive(context, subgraph, complement, excluded_ext,
									complement_neighborhood);
	}
}

/*
 * Expand 'subgraph' using it's neighborhood and try to find complement for it
 */
static void
enumerate_csg_recursive(DPHypContext *context, HyperNode *subgraph,
						bitmapword excluded, bitmapword subgraph_neighborhood)
{
	SubsetIteratorState subset_iter;

	Assert(!bmw_is_empty(subgraph_neighborhood));

	subset_iterator_init(&subset_iter, subgraph_neighborhood);
	while (subset_iterator_next(&subset_iter))
	{
		HyperNode *expanded_subgraph;

		expanded_subgraph = get_hypernode(context, subgraph->set | subset_iter.subset);
		if (hypernode_has_rel(expanded_subgraph))
		{
			bitmapword subgraph_excluded = expanded_subgraph->set | bmw_all_bit_set(bmw_rightmost_one_pos(expanded_subgraph->set));
			bitmapword subnode_neighbors = get_neighbors(context, expanded_subgraph, subgraph_excluded);
			if (!bmw_is_empty(subnode_neighbors))
				emit_csg(context, expanded_subgraph, subnode_neighbors);
		}
	}

	excluded |= subgraph_neighborhood;

	subset_iterator_init(&subset_iter, subgraph_neighborhood);
	while (subset_iterator_next(&subset_iter))
	{
		bitmapword current_neighborhood;
		HyperNode *expanded_subgraph;

		current_neighborhood = get_neighbors_iter(context, subgraph->set | subset_iter.subset,
												  excluded, &subset_iter);
		if (bmw_is_empty(current_neighborhood))
			continue;
		expanded_subgraph = get_hypernode(context, subgraph->set | subset_iter.subset);
		enumerate_csg_recursive(context, expanded_subgraph, excluded,
								current_neighborhood);
	}
}

/* Entry point of DPHyp join search */
static void
solve(DPHypContext *context)
{
	int base_hypernodes_count = list_length(context->base_hypernodes);

	/*
	 * For initial nodes we must iterate backwards to prevent exploring duplicates
	 */
	for (int i = base_hypernodes_count - 1; i >= 0; i--)
	{
		bitmapword neighborhood;
		bitmapword excluded;
		HyperNode *subgraph = (HyperNode *) list_nth(context->base_hypernodes, i);

		excluded = bmw_all_bit_set(i);
		neighborhood = get_neighbors(context, subgraph, excluded);

		if (!bmw_is_empty(neighborhood))
		{
			emit_csg(context, subgraph, neighborhood);
			enumerate_csg_recursive(context, subgraph, excluded, neighborhood);
		}

		/*
		 * Add this in case planning will take too long and user
		 * request cancellation.
		 */
		CHECK_FOR_INTERRUPTS();
	}
}

/*
 * Map Relids specified in 'original' to internal presentation based on id of relation
 */
static bitmapword
map_to_internal_bms(List *initial_rels, Bitmapset *original)
{
	bitmapword target;
	ListCell *lc;
	int i;

	/* 
	 * We must iterate over initial_rels, because RelOptInfo
	 * can represent join, thus it's 'relids' is not singleton.
	 */
	target = 0;
	i = 0;
	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		if (bms_is_subset(rel->relids, original))
			target = bmw_add_member(target, i);
		++i;
	}

	return target;
}

static HyperNode *
get_hypernode(DPHypContext *context, bitmapword nodes)
{
	HyperNode *node;
	bitmapword key = nodes;
	bool found;

	node = hash_search(context->dptable, &key, HASH_ENTER, &found);

	if (!found)
	{
		int idx;

		node->set = nodes;
		node->representative = bmw_first(nodes);
		node->rel = NULL;
		node->candidates = NIL;

		node->simple_edges = 0;
		idx = -1;
		while ((idx = bmw_next_member(nodes, idx)) >= 0)
		{
			node->simple_edges |= context->simple_edges[idx];
		}
		node->simple_edges = bmw_difference(node->simple_edges, node->set);
	}

	return node;
}

/*
 * Get 'RelOptInfo' for given 'HyperNode' and possibly building it.
 * This is called at the end of DPhyp when we are building plan.
 */
static RelOptInfo *
hypernode_get_rel(DPHypContext *context, HyperNode *node)
{
	ListCell *lc;
	RelOptInfo *final_rel;
	HyperNode *left_node;
	HyperNode *right_node;

	/* 
	 * To build final relation we act like vanilla PostgreSQL, but in
	 * top-down approach - build final relation by recursively building
	 * all candidates.
	 */
	if (node->rel != NULL)
		return node->rel;

	final_rel = NULL;
	left_node = NULL;
	right_node = NULL;
	foreach(lc, node->candidates)
	{
		RelOptInfo *left_rel;
		RelOptInfo *right_rel;
		RelOptInfo *join_rel;

		if (left_node == NULL)
		{
			left_node = (HyperNode *)lfirst(lc);
			continue;
		}

		right_node = (HyperNode *)lfirst(lc);

		left_rel = hypernode_get_rel(context, left_node);
		if (left_rel == NULL)
			goto loop_end;

		right_rel = hypernode_get_rel(context, right_node);
		if (right_rel == NULL)
			goto loop_end;

		join_rel = make_join_rel(context->root, left_rel, right_rel);
		if (join_rel == NULL)
			goto loop_end;

		if (final_rel == NULL)
			final_rel = join_rel;

loop_end:
		left_node = NULL;
		right_node = NULL;
	}

	if (final_rel == NULL)
	{
		/*
		 * If we are here, then we are unable to create rel from this node,
		 * then mark this node as invalid to prevent multiple recursive calls
		 * by clearing candidate List.
		 */
		node->candidates = NIL;
		return NULL;
	}

	generate_partitionwise_join_paths(context->root, final_rel);
#if PG_MAJORVERSION_NUM < 16
    if (!bms_equal(context->root->all_baserels, final_rel->relids))
#else	
    if (!bms_equal(context->root->all_query_rels, final_rel->relids))
#endif
		generate_useful_gather_paths(context->root, final_rel, false);
	set_cheapest(final_rel);
	node->rel = final_rel;
	return final_rel;
}

static void
initialize_hypernodes(DPHypContext *context)
{
	ListCell *lc;
	HTAB *dptable;
	HASHCTL hctl;
	int i;

	/* Initial size of HTAB given from 'build_join_rel_hash' */
	hctl.keysize = sizeof(bitmapword);
	hctl.entrysize = sizeof(HyperNode);
	hctl.hash = bmw_hash;
	hctl.match = bmw_match;
	hctl.hcxt = CurrentMemoryContext;
	dptable = (HTAB *)hash_create("DPhyp hypernode table", 256L, &hctl,
								  HASH_ELEM | HASH_FUNCTION | HASH_COMPARE | HASH_CONTEXT);
	i = 0;
	foreach (lc, context->initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		HyperNode *entry;
		bitmapword set;
		bool found;

		set = bmw_make_singleton(i);
		entry = (HyperNode *) hash_search(dptable, &set, HASH_ENTER, &found);

		Assert(!found);

		entry->rel = rel;
		entry->candidates = NIL;
		entry->representative = i;
		entry->set = set;
		entry->simple_edges = context->simple_edges[i];

		context->base_hypernodes = lappend(context->base_hypernodes, entry);
		++i;
	}

	context->dptable = dptable;
}

/* Structure that stores information of Union/Set algorithm */
typedef struct us_state
{
	/* Array of leaders */
	int *leaders;
	/* Array of ranks for each node */
	int *ranks;
	/* Size of 'leaders' and 'ranks' arrays */
	int size;
} us_state;

static void
us_init(us_state *state, int size)
{
	int *leaders;
	int *ranks;

	leaders = palloc(sizeof(int) * size);
	for (size_t i = 0; i < size; i++)
	{
		leaders[i] = i;
	}
	ranks = palloc0(sizeof(int) * size)	;

	state->leaders = leaders;
	state->ranks = ranks;
	state->size = size;
}

static int
us_leader(us_state *state, int node)
{
	Assert(node < state->size);
	if (state->leaders[node] == node)
		return node;
	else
		return state->leaders[node] = us_leader(state, state->leaders[node]);
}

static void
us_union(us_state *state, int a, int b)
{
	int a_leader;
	int b_leader;

	a_leader = us_leader(state, a);
	b_leader = us_leader(state, b);

	if (state->ranks[a_leader] == state->ranks[b_leader])
		state->ranks[a_leader]++;

	if (state->ranks[a_leader] < state->ranks[b_leader])
		state->leaders[a_leader] = b_leader;
	else
		state->leaders[b_leader] = a_leader;
}

static bool us_all_connected(us_state *state)
{
	int prev_leader = -1;
	for (size_t i = 0; i < state->size; i++)
	{
		int leader = us_leader(state, i);
		if (prev_leader == -1)
		{
			prev_leader = leader;
		}
		else if (prev_leader != leader)
		{
			return false;
		}
	}

	return true;
}

static bitmapword *
us_collect(us_state *state, int *size)
{
	bitmapword *disjoint_sets;
	bitmapword *result;
	int result_size;
	int idx;

	disjoint_sets = palloc0(sizeof(bitmapword) * state->size);
	result_size = 0;
	for (size_t i = 0; i < state->size; i++)
	{
		int leader = us_leader(state, i);
		if (bmw_is_empty(disjoint_sets[leader]))
		{
			disjoint_sets[leader] = bmw_make_singleton(i);
			result_size++;
		}
		else
			disjoint_sets[leader] = bmw_add_member(disjoint_sets[leader], i);
	}

	result = palloc(sizeof(bitmapword) * result_size);
	idx = 0;
	for (size_t i = 0; i < state->size; i++)
	{
		if (bmw_is_empty(disjoint_sets[i]))
			continue;

		result[idx] = disjoint_sets[i];
		idx++;
	}

	pfree(disjoint_sets);
	*size = result_size;
	return result;
}

static void
us_free(us_state *state)
{
	pfree(state->leaders);
	pfree(state->ranks);
}

static void
hyperedge_array_add(EdgeArray *array, HyperEdge edge)
{
	int low;
	int high;
	int mid;

	if (array->size == 0)
	{
		/* If array is empty just do allocation and insert edge */
		array->capacity = 4;
		array->size = 1;
		array->edges = palloc(sizeof(HyperEdge) * array->capacity);
		array->edges[0] = edge;
		return;
	}

	/* 
	 * Edge array stored in sorted way for 2 reasons:
	 * 
	 * 1. make use of 'start_index'
	 * 2. prevent duplicates in different hyperedges representation
	 * 
	 * So, here we just perform binary search to find insertion place
	 * and also comparison tells us if there is any duplicate.
	 */
	low = 0;
	high = array->size;
	while (low < high)
	{
		int cmp;
		mid = low + ((high - low) / 2);

		cmp = hyperedge_cmp(edge, array->edges[mid]);

		/* Found duplicate */
		if (cmp == 0)
			return;

		if (cmp < 0)
			high = mid;
		else
			low = mid + 1;
	}

	/* Maybe last element is duplicate */
	if (hyperedge_cmp(edge, array->edges[low]) == 0)
		return;

	/* Suitable position found - adjust edges and insert */
	if (array->size == array->capacity)
	{
		array->capacity *= 2;
		array->edges = repalloc(array->edges, sizeof(HyperEdge) * array->capacity);
	}

	Assert(low <= array->size);
	if (low == array->size)
	{
		array->edges[array->size] = edge;
	}
	else
	{
		memmove(&array->edges[low + 1], &array->edges[low],
				sizeof(HyperEdge) * (array->size - low));
		array->edges[low] = edge;
	}

	array->size++;
}

static inline HyperEdge
hyperedge_swap(HyperEdge edge)
{
	HyperEdge new_edge;
	new_edge.left = edge.right;
	new_edge.right = edge.left;
	return new_edge;
}

static void
distribute_simple_hyperedge(DPHypContext *context, HyperEdge edge)
{
	bitmapword left_bmw;
	bitmapword right_bmw;
	int left_idx;
	int right_idx;

	Assert(hyperedge_is_valid(edge) && hyperedge_is_simple(edge));
	
	left_idx = bmw_first(edge.left);
	right_idx = bmw_first(edge.right);

	left_bmw = context->simple_edges[left_idx];
	right_bmw = context->simple_edges[right_idx];
	context->simple_edges[left_idx] = bmw_add_member(left_bmw, right_idx);
	context->simple_edges[right_idx] = bmw_add_member(right_bmw, left_idx);
}

static void
distribute_hyperedge(DPHypContext *context, HyperEdge edge)
{
	Assert(hyperedge_is_valid(edge));

	if (hyperedge_is_simple(edge))
	{
		distribute_simple_hyperedge(context, edge);
	}
	else
	{
		/* Add hyperedge only to it's representer, not every node in vertexes */
		hyperedge_array_add(&context->complex_edges[bmw_first(edge.left)], edge);
		distribute_cjs(context, edge.left);
		edge = hyperedge_swap(edge);
		hyperedge_array_add(&context->complex_edges[bmw_first(edge.left)], edge);
		distribute_cjs(context, edge.left);
	}
}

static void
distribute_cjs(DPHypContext *context, bitmapword cjs)
{
	HyperEdge edge;
	int idx1;
	int idx2;

	if (bmw_is_empty(cjs) || bmw_is_singleton(cjs))
		return;

	idx1 = -1;
	while ((idx1 = bmw_next_member(cjs, idx1)) >= 0)
	{
		edge.left = bmw_make_singleton(idx1);
		idx2 = idx1;

		while ((idx2 = bmw_next_member(cjs, idx2)) >= 0)
		{
			edge.right = bmw_make_singleton(idx2);
			distribute_simple_hyperedge(context, edge);
		}
	}
}

static bitmapword *
collect_disjoint_sets(DPHypContext *context, int *out_size)
{
	us_state state;
	bitmapword *disjoint_sets;

	us_init(&state, context->edges_size);
	for (int i = 0; i < context->edges_size; ++i)
	{
		bitmapword simple_edge = context->simple_edges[i];
		int idx;

		idx = -1;
		while ((idx = bmw_next_member(simple_edge, idx)) >= 0)
		{
			us_union(&state, i, idx);
		}
	}

	/* As simple heuristic we may find that 'simple_edges' can detect
	 * that all nodes are connected to each other and we can stop now.
	 */
	if (us_all_connected(&state))
	{
		us_free(&state);
		return NULL;
	}

	/*
	 * Disjoint sets exist and be have to generate hyperedges
	 * covering all such disjoint sets.  So process complex edges,
	 * collect disjoint sets and generate hyperedges.
	 */
	for (int i = 0; i < context->edges_size; ++i)
	{
		EdgeArray *edges = &context->complex_edges[i];

		for (int j = 0; j < edges->size; ++j)
		{
			HyperEdge edge = edges->edges[j];
			List *left_vertices;
			int idx;

			idx = -1;
			left_vertices = NIL;
			while ((idx = bmw_next_member(edge.left, idx)) >= 0)
			{
				left_vertices = lappend_int(left_vertices, idx);
			}

			idx = -1;
			while ((idx = bmw_next_member(edge.right, idx)) >= 0)
			{
				ListCell *lc;
				foreach(lc, left_vertices)
				{
					int left_vertex = lfirst_int(lc);
					us_union(&state, left_vertex, idx);
				}
			}
		}
	}

	disjoint_sets = us_collect(&state, out_size);
	us_free(&state);
	if (*out_size <= 1)
	{
		/* All nodes are connected to each other */
		if (disjoint_sets != NULL)
			pfree(disjoint_sets);
		return NULL;
	}

	return disjoint_sets;
}

static List *
collect_disjoint_rels(DPHypContext *context)
{
	List *result;
	int disjoint_sets_size;
	bitmapword *disjoint_sets;

	disjoint_sets = collect_disjoint_sets(context, &disjoint_sets_size);
	if (disjoint_sets == NULL)
		return NIL;

	/* For each disjoint set collect it's RelOptInfo (build lazy) */
	result = NIL;
	for (int i = 0; i < disjoint_sets_size; ++i)
	{
		bitmapword set = disjoint_sets[i];
		RelOptInfo *rel;
		HyperNode *node;

		node = hash_search(context->dptable, &set, HASH_FIND, NULL);
		if (!(node && hypernode_has_rel(node)))
		{
			list_free(result);
			result = NIL;
			break;
		}

		rel = hypernode_get_rel(context, node);
		if (!rel)
		{
			/* This relation is unable to build */
			list_free(result);
			result = NIL;
			break;
		}

		result = lappend(result, rel);
	}

	pfree(disjoint_sets);
	return result;
}

/* 
 * Get start index to iterate over complex edges.
 * Suitable index found using 'excluded' set - first 0 from start.
 */
static int
get_start_index(EdgeArray *edges, bitmapword excluded)
{
	int index;
	int lowest_bit;

	if (edges->start_idx_size == 0)
		return edges->size;

	/*
	 * 'start_idx' primarily used to effectively truncate edges that will not
	 * satisfy 'bmw_overlaps' with 'excluded' set of nodes.
	 * The main observation is that often we have all leading 1 in 'excluded',
	 * so right vertex in any edge with first bit in that range definitely will
	 * return 'false'.
	 * To address this 'start_idx' is used. It is an array:
	 *
	 * [number of leading 0] -> index in 'edges' array
	 *
	 * 'edges' array is sorted by number of leading 0, so we can assert that
	 * if we have 0010 then 0100 will also not overlap with 0001.
	 *
	 * To search suitable position we find first 0 bit after some leading 1.
	 * This is done by inverse - add 1 to sequence of leading 1 and count
	 * produced amount of 0. e.g.
	 *
	 * 1001111 + 1 -> 1010000 (4 leading 1s == 4 leading 0s)
	 * 
	 * This is also used when finding connection between hypernodes.
	 * Difference is that caller must decrement 1 from right side (for
	 * which to check for subset), e.g.
	 * 
	 * 01010100 - 1 -> 01010011
	 */
	Assert(excluded != ~((bitmapword)0));
	lowest_bit = bmw_rightmost_one_pos(excluded + 1);

	if (edges->start_idx_size <= lowest_bit)
		return edges->size;

	index = (int)edges->start_idx[lowest_bit];
	Assert(0 <= index && index < BITS_PER_BITMAPWORD);
	return index;
}

static void
compute_start_index(DPHypContext *context)
{
	for (size_t i = 0; i < context->edges_size; i++)
	{
		EdgeArray *edges = &context->complex_edges[i];
		char prev_idx;
		int prev_lowest;

		if (edges->size == 0)
		{
			edges->start_idx = NULL;
			edges->start_idx_size = 0;
			continue;
		}

		/*
		 * Array indexed by number of bits, so there 2 observations:
		 *
		 * 1. Maximum useful size of this index does not exceed largest
		 *    number of leading bits, so we allocate that amount.
		 *    Array is sorted, so just get size of last hyperedge.
		 * 2. We should reserve special value for 0 number of set bits. This
		 *    value always is 0 (have to traverse all array).
		 */
		edges->start_idx_size = bmw_first(edges->edges[edges->size - 1].right) + 1;
		edges->start_idx = palloc(sizeof(int8) * edges->start_idx_size);

		if (edges->size == 1)
		{
			/*
			 * In case of simple query there may be single complex edge.
			 * You can observe, that this will be array of 0.
			 */
			memset(edges->start_idx, 0, sizeof(int8) * edges->start_idx_size);
			continue;
		}

		/* Set -1 as indicator, that we do not have value set yet */
		memset(edges->start_idx, -1, sizeof(int8) * edges->start_idx_size);

		edges->start_idx[0] = 0;
		prev_lowest = 0;

		/*
		 * Proceed in 2 runs:
		 *
		 * 1. Iterate over all edges and for each possible leading zero bit
		 *    count save position where it starts. Here we use knowledge,
		 *    that hyperedges are sorted, so just track previous 'lowest'
		 *    number and compare with current
		 * 2. Iterate over 'start_index' array and fill missing indexes.
		 *    If value is absent (-1), then set it to previous value (we
		 *    iterate left->right).
		 */

		/* First run - set all possible values */
		for (size_t j = 0; j < edges->size; j++)
		{
			int cur_lowest = bmw_first(edges->edges[j].right);
			if (cur_lowest == prev_lowest)
				continue;

			prev_lowest = cur_lowest;
			edges->start_idx[cur_lowest] = j;
		}

		/* Second run - fill missing indexes */
		prev_idx = 0;
		for (size_t j = 0; j < edges->start_idx_size; j++)
		{
			if (edges->start_idx[j] == -1)
				edges->start_idx[j] = prev_idx;
			else
				prev_idx = edges->start_idx[j];
		}
	}
}

static void
process_edge_bms_pair(DPHypContext *context, Bitmapset *lhs, Bitmapset *rhs)
{
	HyperEdge edge;

	edge.left = map_to_internal_bms(context->initial_rels, lhs);
	if (bmw_is_empty(edge.left))
		return;
	edge.right = map_to_internal_bms(context->initial_rels, rhs);
	if (bmw_is_empty(edge.right))
		return;

	distribute_hyperedge(context, edge);
}

static void
initialize_edges(PlannerInfo *root, List *initial_rels, DPHypContext *context)
{
	ListCell *lc1;
	ListCell *lc2;
	bool has_eclass_joins;

	/* 
	 * Edges for algorithm taken from 3 places:
	 * 
	 * 1. RelOptInfo->joininfo - generic expressions
	 * 2. PlannerInfo->eq_classes - INNER equi-joins
	 * 3. PlannerInfo->join_info_list - non-INNER joins
	 * 
	 * Note that this does not cover all possible edges.
	 * 
	 * Example is 'joininfo' which can have multiple versions for same
	 * expression, but different Relids for left/right side.
	 * This can arise  from non-INNER joins where such redundant RestrictInfo's
	 * store different set of required relations (which turn into vertices of
	 * hyperedges).
	 * When this happens in RestrictInfo->syn_
	 */

	context->edges_size = list_length(initial_rels);
	context->simple_edges = palloc0(sizeof(bitmapword) * list_length(initial_rels));
	context->complex_edges = palloc0(sizeof(EdgeArray) * list_length(initial_rels));

	has_eclass_joins = false;
	foreach(lc1, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc1);

		if (rel->has_eclass_joins)
			has_eclass_joins = true;

		foreach(lc2, rel->joininfo)
		{
			RestrictInfo *rinfo = (RestrictInfo *)lfirst(lc2);

			if (!bms_is_empty(rinfo->left_relids) &&
				!bms_is_empty(rinfo->right_relids) &&
				!bms_overlap(rinfo->left_relids, rinfo->right_relids))
			{
				process_edge_bms_pair(context, rinfo->left_relids, rinfo->right_relids);
			}
			else
			{
				/*
				 * For CJS we must generate all pairs of simple hypernodes
				 */
				bitmapword required_nodes = map_to_internal_bms(initial_rels, rinfo->required_relids);

				if (bmw_is_empty(required_nodes) || bmw_is_singleton(required_nodes))
					continue;

				distribute_cjs(context, required_nodes);
			}
		}
	}

	if (has_eclass_joins)
	{
		/*
		 * Now, we must traverse through all eclasses that can be used as join
		 * clauses and generate edges for them
		 */
		foreach(lc1, root->eq_classes)
		{
			EquivalenceClass *eclass = (EquivalenceClass *)lfirst(lc1);
			bitmapword *eclass_nodes;
			int eclass_nodes_size;

			/* There are definitely no join clauses */
			if (bms_membership(eclass->ec_relids) != BMS_MULTIPLE)
				continue;

			eclass_nodes = palloc0(sizeof(bitmapword) * list_length(eclass->ec_members));
			eclass_nodes_size = 0;

			foreach(lc2, eclass->ec_members)
			{
				EquivalenceMember *member = (EquivalenceMember *)lfirst(lc2);

				if (member->em_is_const || bms_is_empty(member->em_relids))
					continue;

				eclass_nodes[eclass_nodes_size] = map_to_internal_bms(initial_rels, member->em_relids);
				if (bmw_is_empty(eclass_nodes[eclass_nodes_size]))
					continue;
				eclass_nodes_size++;
			}

			if (eclass_nodes_size == 0)
			{
				pfree(eclass_nodes);
				continue;
			}

			for (int i = 0; i < eclass_nodes_size; i++)
			{
				bitmapword left = eclass_nodes[i];
				for (int j = i + 1; j < eclass_nodes_size; j++)
				{
					bitmapword right = eclass_nodes[j];

					if (bmw_overlap(left, right))
					{
						distribute_cjs(context, left | right);
					}
					else
					{
						HyperEdge edge;
						edge.left = left;
						edge.right = right;
						distribute_hyperedge(context, edge);
					}
				}
			}

			pfree(eclass_nodes);
		}
	}

	/*
	 * Join order restrictions also impose restrictions on join order.
	 */
	foreach(lc1, root->join_info_list)
	{
		SpecialJoinInfo *sjinfo = (SpecialJoinInfo *)lfirst(lc1);

		process_edge_bms_pair(context, sjinfo->syn_lefthand, sjinfo->syn_righthand);
		process_edge_bms_pair(context, sjinfo->min_lefthand, sjinfo->min_righthand);
	}

	if (dphyp_cj_strategy == CJ_STRATEGY_DETECT)
	{
		/* Generate all possible pairs for each disjoint set */
		bitmapword *disjoint_sets;
		int disjoint_sets_size;

		disjoint_sets = collect_disjoint_sets(context, &disjoint_sets_size);
		if (disjoint_sets != NULL && 1 < disjoint_sets_size)
		{
			for (size_t i = 0; i < disjoint_sets_size - 1; i++)
			{
				bitmapword left = disjoint_sets[i];
				for (size_t j = i + 1; j < disjoint_sets_size; j++)
				{
					bitmapword right = disjoint_sets[j];
					HyperEdge edge;

					edge.left = left;
					edge.right = right;

					distribute_hyperedge(context, edge);
				}
			}
			
			pfree(disjoint_sets);
		}
	}

	compute_start_index(context);
}

static RelOptInfo *
dphyp(DPHypContext *context, PlannerInfo *root, List *initial_rels)
{
	HyperNode *result;
	bitmapword all_query_nodes;

	context->initial_rels = initial_rels;
	context->root = root;
	context->base_hypernodes = NIL;
	initialize_edges(root, initial_rels, context);

	initialize_hypernodes(context);

	solve(context);

	all_query_nodes = bmw_all_bit_set(list_length(initial_rels) - 1);
	result = hash_search(context->dptable, &all_query_nodes, HASH_FIND, NULL);
	if (!(result && result->candidates))
		return NULL;

	return hypernode_get_rel(context, result);
}

static RelOptInfo *
dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	RelOptInfo *rel;
	List *saved_join_rel_list;
	DPHypContext context;
	List *disjoint_rels;

	if (!dphyp_enabled ||
		BITS_PER_BITMAPWORD <= list_length(initial_rels))
	{
		if (prev_join_search_hook)
			return prev_join_search_hook(root, levels_needed, initial_rels);
		if (enable_geqo && levels_needed >= geqo_threshold)
			return geqo(root, levels_needed, initial_rels);
		return standard_join_search(root, levels_needed, initial_rels);
	}

	/*
	 * Before proceeding we store 'join_rel_list' in case we fail
	 * to find query plan.  Restore it before proceeding to DPsize otherwise
	 * it will throw "failed to build %d-way join"
	 */
	saved_join_rel_list = list_copy(root->join_rel_list);

	rel = dphyp(&context, root, initial_rels);
	statistics.total_runs++;

	/* Successfully found join order */
	if (rel)
	{
		list_free(saved_join_rel_list);
		return rel;
	}

	statistics.failed_runs++;

	if ((dphyp_cj_strategy == CJ_STRATEGY_PASS && (disjoint_rels = collect_disjoint_rels(&context)) != NIL))
	{
		initial_rels = disjoint_rels;
		levels_needed = list_length(disjoint_rels);
	}
	else
	{
		/* Restore state before proceeding to DPsize/GEQO */
		list_free(root->join_rel_list);
		root->join_rel_list = saved_join_rel_list;
	}

	/* Fallback to conventional DPsize/GEQO */
	if (prev_join_search_hook)
		return prev_join_search_hook(root, levels_needed, initial_rels);
	if (enable_geqo && levels_needed >= geqo_threshold)
		return geqo(root, levels_needed, initial_rels);
	return standard_join_search(root, levels_needed, initial_rels);
}

void
_PG_init(void)
{
	DefineCustomBoolVariable("pg_dphyp.enabled",
							 "pg_dphyp join enumeration algorithm is enabled",
							 NULL,
							 &dphyp_enabled,
							 dphyp_enabled,
							 PGC_USERSET,
							 0, NULL, NULL, NULL);
	DefineCustomEnumVariable("pg_dphyp.cj_strategy",
							 "Specifies how extension handles disjoint relations "
							 "that arise from CROSS JOINs",
							 NULL,
							 &dphyp_cj_strategy,
							 dphyp_cj_strategy,
							 cross_join_strategy_options,
							 PGC_USERSET,
					 		 0, NULL, NULL, NULL);

#if PG_MAJORVERSION_NUM < 15
	EmitWarningsOnPlaceholders("pg_dphyp");
#else
	MarkGUCPrefixReserved("pg_dphyp");
#endif

	prev_join_search_hook = join_search_hook;
	join_search_hook = dphyp_join_search;
}

void
_PG_fini(void)
{
	join_search_hook = prev_join_search_hook;
}
