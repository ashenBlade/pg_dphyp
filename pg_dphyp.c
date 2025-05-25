#include "postgres.h"
#include "fmgr.h"
#include "nodes/bitmapset.h"
#include "optimizer/geqo.h"
#include "optimizer/paths.h"
#include "optimizer/pathnode.h"
#include "utils/builtins.h"
#include "utils/hsearch.h"
#include "utils/guc.h"
#include "miscadmin.h"
#include "limits.h"

#include "simplebms.h"

PG_MODULE_MAGIC;

static join_search_hook_type prev_join_search_hook = NULL;

/* GUC */
/* Extension is enabled and should run DPhyp */
static bool dphyp_enabled = true;
/* Do not apply DPhyp if SJ found (LEFT/RIGHT/OUTER etc...) */
static bool dphyp_skip_sj = false;

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
	bitmapword nodes;

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

typedef struct HyperEdge
{
	bitmapword left;
	bitmapword right;
} HyperEdge;

typedef struct EdgeArray
{
	int capacity;
	int size;
	HyperEdge *edges;
} EdgeArray;

static inline bool hyperedge_is_simple(HyperEdge edge)
{
	return !bmw_is_singleton(bmw_union(edge.left, edge.right));
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

static inline int hyperedge_cmp(HyperEdge a, HyperEdge b)
{
	/* 
	 * Simple integer tuple comparison
	 */
	bitmapword t = a.left - b.left;
	if (t != 0)
		return t;
	t = a.right - b.right;
	return t;
}

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

	/* 
	 * Cached bitmap of all nodes appearing in query.
	 * Used in 'emit_csg_cmp' when checking rel before
	 * calling 'generate_useful_gather_paths'.
	 */
	bitmapword all_query_nodes;
} DPHypContext;


static HTAB *create_dptable(List *base_hypernodes);
static HyperNode *create_initial_hypernode(PlannerInfo *root, RelOptInfo *rel,
										   int id, DPHypContext *context);
static bitmapword map_to_internal_bms(List *initial_rels, Bitmapset *original);
static void initialize_edges(PlannerInfo *root, List *initial_rels,
							 DPHypContext *context);
static void distribute_cjs(DPHypContext *context, bitmapword cjs);
static void distribute_hyperedge(DPHypContext *context, HyperEdge edge);
static HyperEdge hyperedge_swap(HyperEdge edge);
static void hyperedge_array_add(EdgeArray *array, HyperEdge edge);
static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes);
static bitmapword get_neighbors(DPHypContext *context, HyperNode *node, bitmapword excluded);
static bool hypernode_has_direct_edge_with(DPHypContext *context, HyperNode *node, int id);
static bool hypernode_has_edge_with(DPHypContext *context, HyperNode *node, bitmapword bms);
static void emit_csg_cmp(DPHypContext *context, HyperNode *subgroup,
						 HyperNode *complement);
static void enumerate_cmp_recursive(DPHypContext *context, HyperNode *node,
									HyperNode *complement, bitmapword excluded);
static void emit_csg(DPHypContext *context, HyperNode *node);
static void enumerate_csg_recursive(DPHypContext *context, HyperNode *node,
									bitmapword excluded);
static void solve(DPHypContext *context);
static RelOptInfo *dphyp(PlannerInfo *root, List *initial_rels);
static RelOptInfo *hypernode_get_rel(DPHypContext *context, HyperNode *node);
static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed,
									 List *initial_rels);

/* 
 * Structure used as state for enumerating subsets of given bitmap
 */
typedef struct
{
	/*
	 * Current subset to return. 0 means no more subsets.
	 */
	bitmapword state;
	/* 
	 * Initial bitmap that used as mask to iterate.
	 */
	bitmapword init;
} SubsetIteratorState;
static void subset_iterator_init(SubsetIteratorState *state, bitmapword bmw);
static bool subset_iterator_next(SubsetIteratorState *state, bitmapword *result);

/* 
 * Check that we calculated any query plan for this hypernode
 */
static inline bool hypernode_has_rel(HyperNode *node)
{
	return node->rel != NULL || node->candidates != NIL;
}

/* 
 * Check that query contains any special joins (LEFT/ANTI/SEMI etc...)
 */
static inline bool contains_sj(PlannerInfo *root)
{
	/* 
	 * There may be more fine-grained checks, but for small OLTP queries
	 * this will introduce too much overhead.
	 * So, just check whether or not we have SJ in query at all.
	 */
	return root->join_info_list != NIL;
}


/*
 * Get gitmap of neighbors for node excluding all specified.
 * Corresponds to 'N(S, X)' function in paper.
 */
static bitmapword
get_neighbors(DPHypContext *context, HyperNode *node, bitmapword excluded)
{
	bitmapword neighbors;
	int idx;
	
	/* Collect neighbors from simple edges */
	neighbors = node->simple_edges;

	/* And then from complex only if they do not overlap with excluded */
	idx = -1;
	while ((idx = bmw_next_member(node->nodes, idx)) >= 0)
	{
		EdgeArray *edges = &context->complex_edges[idx];
		for (int i = 0; i < edges->size; i++)
		{
			HyperEdge edge = edges->edges[i];
			if ( bmw_is_subset(edge.left, node->nodes) &&
				!bmw_overlap(edge.right, neighbors))
			{
				neighbors = bmw_union(neighbors, bmw_first(edge.right));
			}
		}
	}

	neighbors = bmw_difference(neighbors, excluded);
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
	for (int i = 0; i < edges->size; i++)
	{
		HyperEdge edge = edges->edges[i];
		if (edge.left != right_bmw)
			break;

		if (bmw_is_subset(node->nodes, edge.right))
			return true;
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

	Assert(!bmw_overlap(node->nodes, bmw));

	/* Check that we have simple edges that connect to 'bmw' */
	if (bmw_overlap(node->simple_edges, bmw))
		return true;

	/* Now check any complex edge has connection to 'bmw' */
	idx = -1;
	while ((idx = bmw_next_member(node->nodes, idx)) >= 0)
	{
		EdgeArray *edges = &context->complex_edges[idx];
		
		for (int i = 0; i < edges->size; i++)
		{
			HyperEdge edge = edges->edges[i];
			if (bmw_is_subset(edge.left, node->nodes) && 
				bmw_is_subset(edge.right, bmw))
				return true;
		}
	}

	return false;
}

static void
subset_iterator_init(SubsetIteratorState *state, bitmapword bmw)
{
	state->init = bmw;
	state->state = (-bmw) & bmw;
}

static bool
subset_iterator_next(SubsetIteratorState *state, bitmapword *result)
{
	if (state->state == 0)
		return false;

	*result = state->state;
	state->state = (state->state - state->init) & state->init;
	return true;
}

static void
emit_csg_cmp(DPHypContext *context, HyperNode *subgroup, HyperNode *complement)
{
	bitmapword nodes;
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
	nodes = bmw_union(subgroup->nodes, complement->nodes);
	hypernode = get_hypernode(context, nodes);
	hypernode->candidates = lappend(hypernode->candidates, list_make2(subgroup, complement));
}

static void
enumerate_cmp_recursive(DPHypContext *context, HyperNode *node, HyperNode *complement, bitmapword excluded)
{
	bitmapword complement_neighbors;
	SubsetIteratorState subset_iter;
	bitmapword subset;

	complement_neighbors = get_neighbors(context, complement, excluded);
	if (bmw_is_empty(complement_neighbors))
		return;

	subset_iterator_init(&subset_iter, complement_neighbors);
	while (subset_iterator_next(&subset_iter, &subset))
	{
		bitmapword neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bmw_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);

		if (hypernode_has_rel(superset_node) &&
			hypernode_has_edge_with(context, node, neighbor_superset))
			emit_csg_cmp(context, node, superset_node);
	}

	excluded = bmw_union(excluded, complement_neighbors);

	complement_neighbors = get_neighbors(context, complement, excluded);
	if (bmw_is_empty(complement_neighbors))
		return;

	subset_iterator_init(&subset_iter, complement_neighbors);
	while (subset_iterator_next(&subset_iter, &subset))
	{
		bitmapword neighbor_superset;
		HyperNode *superset_node;

		neighbor_superset = bmw_union(complement->nodes, subset);
		superset_node = get_hypernode(context, neighbor_superset);

		enumerate_cmp_recursive(context, node, superset_node, excluded);
	}
}

static void
emit_csg(DPHypContext *context, HyperNode *node)
{
	bitmapword excluded;
	bitmapword neighbors;
	int i;

	excluded = bmw_union(node->nodes, bmw_all_bit_set(node->representative));
	neighbors = get_neighbors(context, node, excluded);
	if (bmw_is_empty(neighbors))
		return;

	i = -1;
	while ((i = bmw_prev_member(neighbors, i)) >= 0)
	{
		HyperNode *complement;
		bitmapword excluded_ext;

		complement = (HyperNode *) list_nth(context->base_hypernodes, i);
		
		/*
		 * Here in original paper we create S = {v} and then check that
		 * edge rhs is subset of S.  But as you can see subset of single element
		 * set is that set itself, so we can make optimized searching
		 * for such edge.
		 */
		if (hypernode_has_direct_edge_with(context, node, i))
			emit_csg_cmp(context, node, complement);

		/*
		 * We are iterating backwards on neighbors, so we have to exclude
		 * all nodes lower current, otherwise, we will get duplicates
		 * and execution time will skyrocket.
		 */
		excluded_ext = bmw_union(excluded, bmw_all_bit_set(i));
		enumerate_cmp_recursive(context, node, complement, excluded_ext);
	}
}

static void
enumerate_csg_recursive(DPHypContext *context, HyperNode *node, bitmapword excluded)
{
	SubsetIteratorState subset_iter;
	bitmapword subset;
	bitmapword neighbors;
	bitmapword excluded_ext;

	neighbors = get_neighbors(context, node, excluded);
	if (bmw_is_empty(neighbors))
		return;

	subset_iterator_init(&subset_iter, neighbors);
	while (subset_iterator_next(&subset_iter, &subset))
	{
		bitmapword superset = bmw_union(node->nodes, subset);
		HyperNode *subnode;

		subnode = get_hypernode(context, superset);
		if (hypernode_has_rel(subnode))
			emit_csg(context, subnode);
	}

	excluded_ext = bmw_union(excluded, neighbors);

	subset_iterator_init(&subset_iter, neighbors);
	while (subset_iterator_next(&subset_iter, &subset))
	{
		bitmapword superset = bmw_union(node->nodes, subset);
		HyperNode *subnode;

		subnode = get_hypernode(context, superset);
		enumerate_csg_recursive(context, subnode, excluded_ext);
	}
}

static void
solve(DPHypContext *context)
{
	int base_hypernodes_count = list_length(context->base_hypernodes);

	/* 
	 * For initial nodes we must iterate backwards to 
	 * prevent exploring duplicates
	 */
	for (int i = base_hypernodes_count - 1; i >= 0; i--)
	{
		HyperNode *node = (HyperNode *) list_nth(context->base_hypernodes, i);
		emit_csg(context, node);
		enumerate_csg_recursive(context, node, bmw_all_bit_set(i));

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
	int id;

	target = 0;
	id = 0;
	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		if (bms_is_subset(rel->relids, original))
		{
			target = bmw_add_member(target, id);
		}

		++id;
	}

	return target;
}

static HyperNode *
create_initial_hypernode(PlannerInfo *root, RelOptInfo *rel, int id,
						 DPHypContext *context)
{
	HyperNode *node;

	node = (HyperNode *)palloc(sizeof(HyperNode));

	node->rel = rel;
	node->candidates = NIL;
	node->representative = id;
	node->nodes = bmw_make_singleton(id);
	node->simple_edges = context->simple_edges[id];

	return node;
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

		node->nodes = nodes;
		node->representative = bmw_first(nodes);
		node->rel = NULL;
		node->candidates = NIL;
		
		node->simple_edges = 0;
		idx = -1;
		while ((idx = bmw_next_member(nodes, idx)) >= 0)
		{
			node->simple_edges = bmw_union(node->simple_edges, context->simple_edges[idx]);
		}
		node->simple_edges = bmw_difference(node->simple_edges, node->nodes);
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
	RelOptInfo *rel;

	if (node->rel != NULL)
		return node->rel;

	rel = NULL;
	foreach(lc, node->candidates)
	{
		HyperNode *node_left;
		HyperNode *node_right;
		RelOptInfo *rel_left;
		RelOptInfo *rel_right;
		RelOptInfo *join_rel;
		List *pair = (List *)lfirst(lc);
		
		Assert(list_length(pair) == 2);

		node_left = (HyperNode *)linitial(pair);
		node_right = (HyperNode *)llast(pair);

		rel_left = hypernode_get_rel(context, node_left);
		if (rel_left == NULL)
			continue;
		rel_right = hypernode_get_rel(context, node_right);
		if (rel_right == NULL)
			continue;

		join_rel = make_join_rel(context->root, rel_left, rel_right);
		if (join_rel == NULL)
			continue;
		
		if (rel == NULL)
			rel = join_rel;
	}

	if (rel == NULL)
	{
		/* 
		 * If we are here, then we are unable to create rel from this node,
		 * then mark this node as invalid to prevent multiple recursive calls
		 * by clearing candidate List.
		 */
		node->candidates = NIL;
		return NULL;
	}

	generate_partitionwise_join_paths(context->root, rel);
    if (!bmw_equal(context->all_query_nodes, node->nodes))
		generate_useful_gather_paths(context->root, rel, false);
	set_cheapest(rel);
	node->rel = rel;
	return rel;
}

static HTAB *
create_dptable(List *base_hypernodes)
{
	ListCell *lc;
	HTAB *dptable;
	HASHCTL hctl;

	/* Initial size of HTAB given from 'build_join_rel_hash' */
	hctl.keysize = sizeof(bitmapword);
	hctl.entrysize = sizeof(HyperNode);
	hctl.hash = bmw_hash;
	hctl.match = bmw_match;
	hctl.hcxt = CurrentMemoryContext;
	dptable = (HTAB *)hash_create("DPhypHyperNodeHashTable", 256L, &hctl,
								  HASH_ELEM | HASH_FUNCTION | HASH_COMPARE | HASH_CONTEXT);

	foreach (lc, base_hypernodes)
	{
		HyperNode *node = (HyperNode *)lfirst(lc);
		HyperNode *entry;
		bool found;

		entry = (HyperNode *) hash_search(dptable, &node->nodes, HASH_ENTER, &found);
		Assert(!found);

		entry->rel = node->rel;
		entry->candidates = node->candidates;
		entry->simple_edges = node->simple_edges;
		entry->representative = node->representative;
		entry->nodes = node->nodes;
	}

	return dptable;
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
		array->capacity = 16;
		array->size = 1;
		array->edges = palloc(sizeof(HyperEdge) * array->capacity);
		array->edges[0] = edge;
		return;
	}

	/* Use binary search to quickly find insert position */
	low = 0;
	high = array->size;
	while (low < high)
	{
		int cmp;
		mid = low + ((high - low) / 2);

		cmp = hyperedge_cmp(edge, array->edges[mid]);
		if (cmp == 0)
			return;

		if (cmp < 0)
			high = mid;
		else
			low = mid + 1;
	}

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
		memmove(&array->edges[array->size], &array->edges[array->size + 1],
				sizeof(HyperEdge) * (array->size - low));
		array->edges[low] = edge;
	}

	array->size++;
}

static HyperEdge
hyperedge_swap(HyperEdge edge)
{
	HyperEdge new_edge;
	new_edge.left = edge.right;
	new_edge.right = edge.left;
	return new_edge;
}

static void
distribute_hyperedge(DPHypContext *context, HyperEdge edge)
{
	Assert(hyperedge_is_valid(edge));

	if (hyperedge_is_simple(edge))
	{
		int left_idx = bmw_next_member(edge.left, -1);
		int right_idx = bmw_next_member(edge.right, -1);
		bitmapword left_bmw;
		bitmapword right_bmw;


		Assert(left_idx >= 0 && right_idx >= 0);

		left_bmw = context->simple_edges[left_idx];
		right_bmw = context->simple_edges[right_idx];
		context->simple_edges[left_idx] = bmw_add_member(left_bmw, right_idx);
		context->simple_edges[right_idx] = bmw_add_member(right_bmw, left_idx);
	}
	else
	{
		int idx;

		idx = -1;
		while ((idx = bmw_next_member(edge.left, idx)) >= 0)
		{
			EdgeArray *array = &context->complex_edges[idx];
			hyperedge_array_add(array, edge);
		}
		
		idx = -1;
		edge = hyperedge_swap(edge);
		while ((idx = bmw_next_member(edge.left, idx)) >= 0)
		{
			EdgeArray *array = &context->complex_edges[idx];
			hyperedge_array_add(array, edge);
		}
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
			distribute_hyperedge(context, edge);
		}
	}
}

static void
initialize_edges(PlannerInfo *root, List *initial_rels, DPHypContext *context)
{
	ListCell *lc1;
	ListCell *lc2;
	bool has_eclass_joins;

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
				/* 
				 * For binary 
				 */
				HyperEdge edge;
				
				edge.left = map_to_internal_bms(initial_rels, rinfo->left_relids);
				if (bmw_is_empty(edge.left))
					continue;
				edge.right = map_to_internal_bms(initial_rels, rinfo->right_relids);
				if (bmw_is_empty(edge.right))
					continue;

				distribute_hyperedge(context, edge);
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

				if (!bmw_is_singleton(left))
					distribute_cjs(context, eclass_nodes[i]);

				for (int j = i + 1; j < eclass_nodes_size; j++)
				{
					bitmapword right = eclass_nodes[j];

					if (bmw_overlap(left, right))
					{
						distribute_cjs(context, bmw_union(left, right));
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
	 * Join order restrictions also impose restrictions on join order
	 */
	foreach(lc1, root->join_info_list)
	{
		SpecialJoinInfo *sjinfo = (SpecialJoinInfo *)lfirst(lc1);
		HyperEdge edge;

		edge.left = map_to_internal_bms(initial_rels, sjinfo->min_lefthand);
		if (bmw_is_empty(edge.left))
			continue;
		edge.right = map_to_internal_bms(initial_rels, sjinfo->min_righthand);
		if (bmw_is_empty(edge.right))
			continue;

		distribute_hyperedge(context, edge);
	}
}

static RelOptInfo *
dphyp(PlannerInfo *root, List *initial_rels)
{
	HTAB *dptable;
	HyperNode *result;
	DPHypContext context;
	List *base_hypernodes;
	ListCell *lc;
	int id;
	bool result_found;

	initialize_edges(root, initial_rels, &context);

	base_hypernodes = NIL;
	id = 0;
	foreach(lc, initial_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		HyperNode *node;

		node = create_initial_hypernode(root, rel, id, &context);
		base_hypernodes = lappend(base_hypernodes, node);
		++id;
	}

	dptable = create_dptable(base_hypernodes);

	context.dptable = dptable;
	context.initial_rels = initial_rels;
	context.root = root;
	context.base_hypernodes = base_hypernodes;
	context.all_query_nodes = bmw_all_bit_set(list_length(initial_rels) - 1);
	solve(&context);
	
	result = hash_search(dptable, &context.all_query_nodes, HASH_FIND, &result_found);
	if (!(result && result->candidates))
		return NULL;

	return hypernode_get_rel(&context, result);
}

static RelOptInfo *
dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	RelOptInfo *rel;
	List *saved_join_rel_list;

	if (!dphyp_enabled ||
		BITS_PER_BITMAPWORD <= levels_needed ||
		(dphyp_skip_sj && contains_sj(root)))
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

	rel = dphyp(root, initial_rels);
	
	/* Successfully found join order */
	if (rel)
	{
		list_free(saved_join_rel_list);
		return rel;
	}

	/* Restore state before proceeding to DPsize/GEQO */
	list_free(root->join_rel_list);
	root->join_rel_list = saved_join_rel_list;

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
	DefineCustomBoolVariable("pg_dphyp.skip_sj",
							 "Do not run DPhyp if any special join detected",
							 NULL,
							 &dphyp_skip_sj,
							 dphyp_skip_sj,
							 PGC_USERSET,
							 0, NULL, NULL, NULL);
	MarkGUCPrefixReserved("pg_dphyp");

	prev_join_search_hook = join_search_hook;
	join_search_hook = dphyp_join_search;
}

void
_PG_fini(void)
{
	join_search_hook = prev_join_search_hook;
}
