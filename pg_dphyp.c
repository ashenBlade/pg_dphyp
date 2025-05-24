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
/* Threshold for GEQO used by DPhyp */
static int dphyp_geqo_threshold = 12;

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
	 * List of hyperedges, that this node belongs to.
	 * Each hypernode is represented as List of nodes and each node is bitmapword - Cross Join Set.
	 * Also, each hyperedge not only contains 2 nodes, but any count.
	 */
	List *hyperedges;
} HyperNode;

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
										   int id, List *initial_rels);
static bitmapword map_to_internal_bms(List *initial_rels, Bitmapset *original);

static int bmw_list_comparator(const ListCell *a, const ListCell *b);
static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes);
static bitmapword get_neighbors(HyperNode *node, bitmapword excluded);
static bool hypernode_has_direct_edge_with(HyperNode *node, int id);
static bool hypernode_has_edge_with(HyperNode *node, bitmapword bms);
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
get_neighbors(HyperNode *node, bitmapword excluded)
{
	/*
	 * Neighbors of hypernode are just all reachable nodes from this node
	 * except excluded.
	 * In such setting, all we need to do is to iterate through all edges
	 * and collect neighbors.
	 */
	ListCell *lc;
	bitmapword neighbors;

	neighbors = 0;
	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;

		foreach (lc2, edge)
		{
			bitmapword set = (bitmapword )lfirst(lc2);

			/*
			 * With concept of Cross Join Sets logic changes a bit - we must
			 * separately handle CJS that contains current node from those that don't.
			 *
			 * If current edge node contains us, then we implicitly join ourself
			 * with every node in that node - treat as multiple CROSS JOIN edges.
			 * So calculate difference '(edge_node / node) / excluded' and add all to our neighbors.
			 *
			 * Otherwise, treat as plain usual edge just like in paper - check
			 * it does not contain any excluded nodes and add only representative.
			 */
			if (bmw_overlap(set, node->nodes))
			{
				bitmapword me_excluded;
				bitmapword without_excluded;

				/* Excluded myself from edge */
				me_excluded = bmw_difference(set, node->nodes);
				if (bmw_is_empty(me_excluded))
					continue;

				/* Exclude excluded nodes from what's left */
				without_excluded = bmw_difference(me_excluded, excluded);
				if (bmw_is_empty(without_excluded))
					continue;

				/* Add rest nodes from this CJS */
				neighbors = bmw_union(neighbors, without_excluded);
			}
			else
			{
				int representative;

				/* Proceed as in original paper - check no excluded and add representative */
				if (bmw_overlap(set, excluded))
					continue;

				Assert(!bmw_is_empty(set));
				representative = bmw_first(set);
				neighbors = bmw_add_member(neighbors, representative);
			}
		}
	}

	return neighbors;
}

/*
 * Check that 'node' has direct edge with node 'id'.
 * This is not the same as 'has_edge_with' because we must check
 * that it has *direct* edge.
 */
static bool
hypernode_has_direct_edge_with(HyperNode *node, int id)
{
	ListCell *lc;
	bitmapword bmw;

	Assert(!bmw_is_member(node->nodes, id));

	bmw = bmw_make_singleton(id);
	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			bitmapword vertex = (bitmapword )lfirst(lc2);

			if (bmw_is_subset(node->nodes, vertex))
			{
				/*
				 * For CJS with us - check 'id' is member of edge
				 */
				if (bmw_is_member(vertex, id))
					return true;
			}
			else
			{
				/*
				 * For usual edges - check this is single element
				 */
				if (bmw_equal(vertex, bmw))
					return true;
			}
		}
	}

	return false;
}

/*
 * Check that 'node' has any edge that can be used as connection to 'bms'
 */
static bool
hypernode_has_edge_with(HyperNode *node, bitmapword bms)
{
	ListCell *lc;

	Assert(!bmw_overlap(node->nodes, bms));

	foreach (lc, node->hyperedges)
	{
		List *edge = (List *)lfirst(lc);
		ListCell *lc2;
		foreach (lc2, edge)
		{
			bitmapword vertex = (bitmapword )lfirst(lc2);

			/*
			 * Of course, we may have vertex intersecting with node,
			 * but logic is the same for both cases - check subset.
			 */
			if (bmw_is_subset(bms, vertex))
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

	complement_neighbors = get_neighbors(complement, excluded);
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
			hypernode_has_edge_with(node, neighbor_superset))
			emit_csg_cmp(context, node, superset_node);
	}

	excluded = bmw_union(excluded, complement_neighbors);

	complement_neighbors = get_neighbors(complement, excluded);
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
	neighbors = get_neighbors(node, excluded);
	if (bmw_is_empty(neighbors))
		return;

	i = -1;
	while ((i = bmw_prev_member(neighbors, i)) >= 0)
	{
		HyperNode *complement;

		/*
		 * Here in original paper we create set with single element.
		 * Then during edges iteration we find edge right node of which is 'subset' of set created above.
		 * But the only possible option - right node is the same set, because subset of single set - it is same same set.
		 */
		complement = (HyperNode *) list_nth(context->base_hypernodes, i);

		if (hypernode_has_direct_edge_with(node, i))
			emit_csg_cmp(context, node, complement);

		/*
		 * We are iterating backwards, so we have to excluded all neighbors that are going to be visited.
		 * Otherwise, we will get duplicates and execution time will skyrocket.
		 */
		enumerate_cmp_recursive(context, node, complement, bmw_union(excluded, bmw_all_bit_set(i)));
	}
}

static void
enumerate_csg_recursive(DPHypContext *context, HyperNode *node, bitmapword excluded)
{
	SubsetIteratorState subset_iter;
	bitmapword subset;
	bitmapword neighbors;
	bitmapword excluded_ext;

	neighbors = get_neighbors(node, excluded);
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

	/* For initial nodes we must iterate backwards to prevent exploring duplicates */
	for (int i = base_hypernodes_count - 1; i >= 0; i--)
	{
		HyperNode *node = (HyperNode *) list_nth(context->base_hypernodes, i);
		emit_csg(context, node);
		enumerate_csg_recursive(context, node, bmw_all_bit_set(i));

		/*
		 * When amount of relation is high, execution time will grow exponentially and user may request cancellation.
		 * This seems a good place to check this.
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
create_initial_hypernode(PlannerInfo *root, RelOptInfo *rel, int id, List *initial_rels)
{
	HyperNode *node;
	List *edges;
	ListCell *lc;

	node = (HyperNode *)palloc(sizeof(HyperNode));

	edges = NIL;

	foreach(lc, rel->joininfo)
	{
		RestrictInfo *rinfo = (RestrictInfo *)lfirst(lc);
		ListCell *lc2;
		bool can_add;
		List *edge = NIL;

		/*
		 * Such check for binary operator is presented in 'make_restrictinfo_internal'
		 */
		if (!bms_is_empty(rinfo->left_relids) &&
			!bms_is_empty(rinfo->right_relids) &&
			!bms_overlap(rinfo->left_relids, rinfo->right_relids))
		{
			/* Create usual edge with 2 relations on both sides */
			bitmapword left_bmw;
			bitmapword right_bmw;

			left_bmw = map_to_internal_bms(initial_rels, rinfo->left_relids);
			if (bmw_is_empty(left_bmw))
				continue;
			right_bmw = map_to_internal_bms(initial_rels, rinfo->right_relids);
			if (bmw_is_empty(right_bmw))
				continue;

			/* For convenience, we will sort vertices based on their representative */
			if (left_bmw < right_bmw)
				edge = list_make2((void *)left_bmw, (void *)right_bmw);
			else
				edge = list_make2((void *)right_bmw, (void *)left_bmw);
		}
		else
		{
			/* This is Cross Join Set with single vertex */
			bitmapword bms;
			bms = map_to_internal_bms(initial_rels, rinfo->required_relids);
			if (!bmw_is_empty(bms))
				edge = list_make1((void *)bms);
		}

		if (!edge)
			continue;

		/*
		 * Verify, that there are no duplicates of edges we are going to add.
		 * This can be true in cases when there are multiple join clauses on same set of relations.
		 * For algorithm purposes, we only need to have one copy of such edge.
		 */
		can_add = true;
		foreach (lc2, edges)
		{
			List *existing_edge = (List *)lfirst(lc2);
			if (list_length(existing_edge) != list_length(edge))
				continue;

			if (list_length(existing_edge) == 1)
			{
				bitmapword existing_vertex = (bitmapword)linitial(existing_edge);
				if (existing_vertex == (bitmapword)linitial(edge))
				{
					can_add = false;
					break;
				}
			}
			else
			{
				bitmapword left_existing = (bitmapword)linitial(existing_edge);
				bitmapword right_existing = (bitmapword)llast(existing_edge);
				bitmapword left_new = (bitmapword)linitial(edge);
				bitmapword right_new = (bitmapword)llast(edge);

				Assert(list_length(existing_edge) == 2);

				/*
				 * Vertices are sorted, so there is only single check is required.
				 */
				if (left_existing == left_new && right_existing == right_new)
				{
					can_add = false;
					break;
				}
			}
		}

		if (can_add)
			edges = lappend(edges, edge);
	}

	/*
	 * PostgreSQL has equivalence classes mechanism that also used for join clauses.
	 * For such cases we can not create conventional edges (1 or 2 vertices) - amount of them will grow exponentially.
	 * So, we extend our hyperedges that way it can contain any amount of vertices.
	 * Each vertex in such hyperedge can be processed same way as regular vertex.
	 */
	if (rel->has_eclass_joins)
	{
		int i = -1;
		while ((i = bms_next_member(rel->eclass_indexes, i)) >= 0)
		{
			EquivalenceClass *ec = (EquivalenceClass *)list_nth(root->eq_classes, i);
			List *edge;

			/* Process only EC containing possible JOIN clauses */
			if (bms_equal(ec->ec_relids, rel->relids))
				continue;

			edge = NIL;
			foreach(lc, ec->ec_members)
			{
				EquivalenceMember *em = (EquivalenceMember *)lfirst(lc);
				bitmapword bmw;

				if (em->em_is_const)
					continue;

				Assert(!bms_is_empty(em->em_relids));
				bmw = map_to_internal_bms(initial_rels, em->em_relids);
				if (!bmw_is_empty(bmw))
					edge = lappend(edge, (void *)bmw);
			}

			if (edge != NIL)
				edges = lappend(edges, edge);
		}
	}

	/* 
	 * OUTER JOINs impose restrictions for join order. Notably, even if
	 * we have some join clause 'a.x = b.x' this does not mean we have to
	 * join only 'a' and 'b', because this can be RIGHT JOIN clause and
	 * then we must join 'b' with 'a' + all relations on left part.
	 * So, all SpecialJoinInfo create hyperedge.
	 */
	foreach (lc, root->join_info_list)
	{
		SpecialJoinInfo *sjinfo = (SpecialJoinInfo *)lfirst(lc);
		bitmapword left_bmw;
		bitmapword right_bmw;
		if (sjinfo->jointype == JOIN_INNER)
			continue;
		
		if (!(bms_is_subset(rel->relids, sjinfo->min_lefthand) ||
			  bms_is_subset(rel->relids, sjinfo->min_righthand)))
			continue;
		
		left_bmw = map_to_internal_bms(initial_rels, sjinfo->min_lefthand);
		if (bmw_is_empty(left_bmw))
			continue;
		right_bmw = map_to_internal_bms(initial_rels, sjinfo->min_righthand);
		if (bmw_is_empty(right_bmw))
			continue;
		edges = lappend(edges, list_make2((void *)left_bmw, (void *)right_bmw));
	}

	node->hyperedges = edges;
	node->rel = rel;
	node->candidates = NIL;
	node->representative = id;
	node->nodes = bmw_make_singleton(id);

	return node;
}

static int
bmw_list_comparator(const ListCell *a, const ListCell *b)
{
	bitmapword a_nodes = (bitmapword)lfirst(a);
	bitmapword b_nodes = (bitmapword)lfirst(b);
	return a_nodes - b_nodes;
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
		List *hyperedges;
		int id;

		node->nodes = nodes;
		node->representative = bmw_first(nodes);
		node->rel = NULL;
		node->candidates = NIL;

		hyperedges = NIL;

		/*
		 * If this hypernode created for the first time, we must initialize it's hyperedges.
		 * This is done by merging all hyperedges from base hypernodes it contains.
		 * When iterating over hyperedges we skip all hypernodes that fully enclosed in current hypernode.
		 *
		 * There are some optimizations for merging base hyperedges.
		 *
		 * 1. We also try to eliminate all duplicates, otherwise the same edge will be created for every base hypernode it mentions.
		 * This is done by sorting all hypernodes in hyperedges and then comparing newly created edge with all added ones.
		 *
		 * 2. All different vertexes fully enclosed in newly created hypernode will be merged into single vertex.
		 * This helps us to reduce vertexes count and thus improve performance.
		 * Note, that this is helpful primarily for EC hyperedges.
		 */
		id = -1;
		while ((id = bmw_next_member(nodes, id)) >= 0)
		{
			HyperNode *base_node = list_nth(context->base_hypernodes, id);
			ListCell *lc;

			foreach (lc, base_node->hyperedges)
			{
				List *edge = (List *)lfirst(lc);
				ListCell *lc2;
				List *new_edge = NIL;
				bitmapword all_my_vertexes = 0;
				bitmapword all_vertexes = 0;
				bool can_add;

				foreach(lc2, edge)
				{
					bitmapword vertex = (bitmapword)lfirst(lc2);

					all_vertexes = bmw_union(all_vertexes, vertex);

					if (bmw_is_subset(vertex, nodes))
						all_my_vertexes = bmw_union(all_my_vertexes, vertex);
					else
						new_edge = lappend(new_edge, (void *) vertex);
				}

				/*
				 * Skip edge that fully enclosed in current hypernode.
				 * It is not useful for further processing, because no neighbors will be detected from it.
				 */
				if (bmw_is_subset(all_vertexes, nodes))
				{
					list_free(new_edge);
					continue;
				}

				if (new_edge == NIL)
					continue;

				/*
				 * Merge all vertexes that belongs to new hypernode into single vertex.
				 */
				if (!bmw_is_empty(all_my_vertexes))
					new_edge = lappend(new_edge, (void *) all_my_vertexes);

				/*
				 * Sort vertexes to further detect duplicates.
				 */
				if (1 < list_length(new_edge))
				{
					/* Ради частного случая в 2 вершины не стоит запускать list_sort */
					if (list_length(new_edge) == 2)
					{
						bitmapword first_vertex = (bitmapword)linitial(new_edge);
						bitmapword second_vertex = (bitmapword)llast(new_edge);
						if (second_vertex < first_vertex)
						{
							linitial(new_edge) = (void *) second_vertex;
							llast(new_edge) = (void *) first_vertex;
						}
					}
					else
					{
						list_sort(new_edge, bmw_list_comparator);
					}
				}

				can_add = true;
				foreach(lc2, hyperedges)
				{
					List *other_edge = (List *)lfirst(lc2);
					bool are_equal;

					/* Edges with different vertexes count definitely are not equal */
					if (list_length(new_edge) != list_length(other_edge))
						continue;

					/* For performance reasons split cases 1, 2 or multiple vertexes */
					are_equal = true;
					if (list_length(new_edge) == 1)
					{
						bitmapword v1 = (bitmapword)linitial(new_edge);
						bitmapword v2 = (bitmapword)linitial(other_edge);
						if (v1 != v2)
							are_equal = false;
					}
					else if (list_length(new_edge) == 2)
					{
						bitmapword first_vertex = (bitmapword)linitial(new_edge);
						bitmapword second_vertex = (bitmapword)llast(new_edge);
						bitmapword other_first_vertex = (bitmapword)linitial(other_edge);
						bitmapword other_second_vertex = (bitmapword)llast(other_edge);
						if (first_vertex != other_first_vertex || second_vertex != other_second_vertex)
							are_equal = false;
					}
					else
					{
						/* Generic case - compare each vertex iteratively */
						ListCell *lc3;
						ListCell *lc4;
						forboth(lc3, new_edge, lc4, other_edge)
						{
							bitmapword first_vertex = (bitmapword)lfirst(lc3);
							bitmapword second_vertex = (bitmapword)lfirst(lc4);
							if (first_vertex != second_vertex)
							{
								are_equal = false;
								break;
							}
						}

					}

					if (are_equal)
					{
						can_add = false;
						break;
					}
				}

				if (can_add)
					hyperedges = lappend(hyperedges, new_edge);
			}
		}

		node->hyperedges = hyperedges;
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
		entry->hyperedges = node->hyperedges;
		entry->representative = node->representative;
		entry->nodes = node->nodes;
	}

	return dptable;
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

	base_hypernodes = NIL;
	id = 0;
	foreach(lc, initial_rels)
	{
		/* 
		 * TODO: create hyperedges and distributes them among hypernodes
		 * 		 instead of creating their own copy for each hypernode.
		 */
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		HyperNode *node;

		node = create_initial_hypernode(root, rel, id, initial_rels);
		base_hypernodes = lappend(base_hypernodes, node);
		++id;
	}

	dptable = create_dptable(base_hypernodes);

	/*
	 * Access paths for all rels in 'initial_rels' already found.
	 * Ready to run DPHyp.
	 */
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
	DefineCustomIntVariable("pg_dphyp.geqo_threshold",
							"Sets the threshold of FROM items beyond which DPhyp algorithm will not be used.",
							NULL,
							&dphyp_geqo_threshold,
							dphyp_geqo_threshold,
							2, BITS_PER_BITMAPWORD - 1,
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
