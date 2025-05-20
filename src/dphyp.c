#include "postgres.h"
#include "utils/hsearch.h"
#include "optimizer/paths.h"
#include "optimizer/pathnode.h"
#include "miscadmin.h"

#include "dphyp.h"
#include "unionset.h"
#include "simplebms.h"

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
	 * May be NULL in cases it is created on demand and used transiently (to pass
	 * as argument)
	 */
	RelOptInfo *rel;

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


	bitmapword all_query_nodes;
} DPHypContext;

static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes);

static HTAB *create_dptable(List *base_hypernodes)
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
		entry->hyperedges = node->hyperedges;
		entry->representative = node->representative;
		entry->nodes = node->nodes;
	}

	return dptable;
}

/**
 * Get gitmap of neighbors for node excluding all specified.
 * Corresponds to 'N(S, X)' function in paper.
 */
static bitmapword get_neighbors(HyperNode *node, bitmapword excluded)
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
static bool hypernode_has_direct_edge_with(HyperNode *node, int id)
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
static bool hypernode_has_edge_with(HyperNode *node, bitmapword bms)
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

typedef struct
{
	bitmapword state;
	bitmapword init;
} SubsetIteratorState;

static void subset_iterator_init(SubsetIteratorState *state, bitmapword bmw)
{
	state->init = bmw;
	state->state = (-bmw) & bmw;
}

static bool subset_iterator_next(SubsetIteratorState *state, bitmapword *result)
{
	if (state->state == 0)
		return false;

	*result = state->state;
	state->state = (state->state - state->init) & state->init;
	return true;
}

static void emit_csg_cmp(DPHypContext *context, HyperNode *subgroup, HyperNode *complement)
{
	bitmapword nodes;
	RelOptInfo *joinrel;
	HyperNode *hypernode;

	/*
	 * For this moment, RelOptInfo for both subgroup and complement must exist
	 */
	Assert(subgroup->rel);
	Assert(complement->rel);

	nodes = bmw_union(subgroup->nodes, complement->nodes);
	hypernode = get_hypernode(context, nodes);

	/*
	 * Our DPhyp logic reuse code from original postgres' DPsize algorithm.
	 * And main function - is 'make_join_rel' which actually creates 'RelOptInfo'
	 * and find all possible paths.
	 *
	 * Current implementation works well with INNER JOIN, but in case of
	 * outer joins (LEFT/SEMI/ANTI/OUTER etc...) it can misbehave due to
	 * join ordering restrictions imposed by outer joints.
	 * In such cases we can not create 'RelOptInfo' and 'make_join_rel' will
	 * just return NULL.
	 */
	joinrel = make_join_rel(context->root, subgroup->rel, complement->rel);
	if (!joinrel)
		return;
		
	/*
	 * Find best path for this rel or update existing one.
	 * Code copied from 'standard_join_search'.
	 */
	generate_partitionwise_join_paths(context->root, joinrel);
    if (!bmw_equal(context->all_query_nodes, hypernode->nodes))
		generate_useful_gather_paths(context->root, joinrel, false);
	set_cheapest(joinrel);

	if (hypernode->rel == NULL)
		hypernode->rel = joinrel;
}

static void enumerate_cmp_recursive(DPHypContext *context, HyperNode *node, HyperNode *complement, bitmapword excluded)
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

		if (superset_node->rel != NULL &&
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

static void emit_csg(DPHypContext *context, HyperNode *node)
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

static void enumerate_csg_recursive(DPHypContext *context, HyperNode *node, bitmapword excluded)
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
		if (subnode->rel != NULL)
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

static void solve(DPHypContext *context)
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
static bitmapword map_to_internal_bms(List *initial_rels, Bitmapset *original)
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

static HyperNode *create_initial_hypernode(PlannerInfo *root, RelOptInfo *rel, int id, List *initial_rels)
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

	node->hyperedges = edges;
	node->rel = rel;
	node->representative = id;
	node->nodes = bmw_make_singleton(id);

	return node;
}

static int bmw_list_comparator(const ListCell *a, const ListCell *b)
{
	bitmapword a_nodes = (bitmapword)lfirst(a);
	bitmapword b_nodes = (bitmapword)lfirst(b);
	return a_nodes - b_nodes;
}

static HyperNode *get_hypernode(DPHypContext *context, bitmapword nodes)
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

static List *collect_disjoint_relations(DPHypContext *context)
{
	/*
	 * When we encounter implicit 'CROSS JOIN's, then our hypergraph will not contain edges between these subgraphs, thus producing disjoint sets of graphs.
	 * DPsize solves this problem just by creating each possible join of relations, but DPhyp will not detect this.
	 *
	 * If there is no 'RelOptInfo' created for top level hypernode, it means we have implicit 'CROSS JOIN's and have to handle it.
	 * We solve this collecting already solved problems (subgraphs) and running DPsize/GEQO with new 'RelOptInfo' set.
	 * To collect all disjoint sets we use UNION-SET structure with ID of hypernode as key for sets.
	 */
	List *disjoint_relations;
	List *disjoint_sets;
	ListCell *lc;
	unionset_state us_state;
	int num_leaders;

	num_leaders = list_length(context->initial_rels);
	us_init(&us_state, num_leaders);

	/* Generate sets */
	foreach(lc, context->base_hypernodes)
	{
		HyperNode *node = (HyperNode *)lfirst(lc);
		ListCell *lc_edge;
		foreach(lc_edge, node->hyperedges)
		{
			ListCell *lc_vertex;
			List *edge = (List *)lfirst(lc_edge);

			foreach(lc_vertex, edge)
			{
				bitmapword vertex = (bitmapword )lfirst(lc_vertex);
				int i = -1;
				while ((i = bmw_next_member(vertex, i)) >= 0)
				{
					if (node->representative != i)
						us_union(&us_state, node->representative, i);
				}
			}
		}
	}

	/* Collect all disjoint sets */
	disjoint_sets = us_collect(&us_state);

	/*
	 * Currently, we can not handle some queries with complex query-tree structures.
	 * i.e. with lots of LEFT/RIGHT/OUTER JOIN's
	 *
	 * For these cases, we return NIL as a signal to run DPsize/GEQO (fallback)
	 */
	Assert(disjoint_sets != NIL);
	if (list_length(disjoint_sets) == 1)
		return NIL;

	/* Iterate through all disjoint sets and collect 'RelOptInfo' for them */
	disjoint_relations = NIL;
	foreach(lc, disjoint_sets)
	{
		HyperNode *hypernode;
		Bitmapset *disjoint_set = (Bitmapset *)lfirst(lc);
		bitmapword simple_disjoint_set = 0;
		int i;

		i = -1;
		while ((i = bms_next_member(disjoint_set, i)) >= 0)
		{
			simple_disjoint_set = bmw_add_member(simple_disjoint_set, i);
		}

		hypernode = get_hypernode(context, simple_disjoint_set);
		if (hypernode->rel == NULL)
		{
			us_free(&us_state);
			list_free(disjoint_sets);
			return NIL;
		}

		disjoint_relations = lappend(disjoint_relations, hypernode->rel);
	}

	list_free(disjoint_sets);
	us_free(&us_state);

	return disjoint_relations;
}


List *dphyp(PlannerInfo *root, int levels_needed, List *initial_rels)
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
	if (result && result->rel)
		return list_make1(result->rel);

	/* If we failed to create plan for whole relation, maybe there are implicit joins */
	return collect_disjoint_relations(&context);
}
