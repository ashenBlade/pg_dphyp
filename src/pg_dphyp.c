#include "postgres.h"
#include "fmgr.h"
#include "nodes/bitmapset.h"
#include "utils/builtins.h"
#include "optimizer/paths.h"
#include "utils/hsearch.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

/* 
 * Represents a hypernode in query graph
 */
typedef struct HyperNode
{
	RelOptInfo *rel;
} HyperNode;

static join_search_hook_type prev_join_search_hook = NULL;

void _PG_init(void);
void _PG_fini(void);

static inline int bms_first(Bitmapset *bms)
{
	int res = bms_next_member(bms, -1);
	Assert(res > 0);
	return res;
}

static Bitmapset *create_bv(RelOptInfo *rel)
{
	Bitmapset *bv;
	int member;

	Assert(bms_membership(rel->relids) == BMS_SINGLETON);
	member = bms_singleton_member(rel->relids);
	bv = bms_make_singleton(member);
	for (int i = 0; i < member; ++i)
	{
		bv = bms_add_member(bv, i);
	}
	
	return bv;
}

static HTAB *create_dptable(List *base_rels)
{
	ListCell *lc;
	HTAB *dptable;
	HASHCTL hctl;

	hctl.keysize = sizeof(Bitmapset *);
	hctl.entrysize = sizeof(RelOptInfo **);
	hctl.hash = bitmap_hash;
	hctl.match = bitmap_match;
	dptable = (HTAB *)hash_create("DPHyp Hash Table", 1024, &hctl,
								  HASH_ELEM | HASH_BLOBS);

	foreach (lc, base_rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		RelOptInfo **entry;
		bool found;

		entry = (RelOptInfo **) hash_search(dptable, rel->relids, HASH_ENTER, &found);
		Assert(!found);

		*entry = rel;
	}

	return dptable;
}

static List *get_neighbours(RelOptInfo *rel, Bitmapset *excluded)
{
	List *neighbours = NIL;
	ListCell *lc;
	foreach(lc, rel->joininfo)
	{
		RestrictInfo *rinfo = (RestrictInfo *)lfirst(lc);
		
	}
}

static void emit_csg_cmp(PlannerInfo *root, RelOptInfo *set, RelOptInfo *complement)
{

}

static void enumerate_csg_cmp_recursive(PlannerInfo *root, RelOptInfo *set, RelOptInfo *complement, Bitmapset *excluded)
{

}

static void emit_csg(PlannerInfo *root, RelOptInfo *set)
{

}

static void enumerate_csg_recursive(PlannerInfo *root, RelOptInfo *set, Bitmapset *excluded)
{

}

static void solve(PlannerInfo *root, List *rels, HTAB *dptable)
{
	RelOptInfo **rels_reversed;
	ListCell *lc;

	rels_reversed = (RelOptInfo **)palloc(list_length(rels) * sizeof(RelOptInfo *));
	
	/* Оборачиваем список, чтобы итерироваться в обратном порядке */
	foreach(lc, rels)
	{
		RelOptInfo *rel = (RelOptInfo *)lfirst(lc);
		rels_reversed[list_length(rels) - foreach_current_index(lc)] = rel;
	}

	for (int i = 0; i < list_length(rels); i++)
	{
		Bitmapset *excluded;
		RelOptInfo *rel = rels_reversed[i];
		emit_csg(root, rel);
		enumerate_csg_recursive(root, rel, create_bv(rel));
	}
}

static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	HTAB *dptable;
	RelOptInfo **entry;
	Bitmapset *bms;

	dptable = create_dptable(initial_rels);

	/* 
	 * Access paths for all rels in 'initial_rels' already found.
	 * Ready to run DPHyp.
	 */

	solve(root, initial_rels, dptable);
	
	elog(ERROR, "DPHyp not implemented");
}

void
_PG_init(void)
{
	/* 
	 * Шаги
	 */
	prev_join_search_hook = join_search_hook;
	join_search_hook = dphyp_join_search;
}

void
_PG_fini(void)
{
	join_search_hook = prev_join_search_hook;
}
