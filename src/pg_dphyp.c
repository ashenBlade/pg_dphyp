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

#include "dphyp.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

static join_search_hook_type prev_join_search_hook = NULL;

/* GUC */
/* Extension is enabled and should run DPhyp */
static bool enabled = true;

void _PG_init(void);
void _PG_fini(void);

static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	List *result;

	if (!enabled || levels_needed <= BITS_PER_BITMAPWORD)
	{
		if (prev_join_search_hook)
			return prev_join_search_hook(root, levels_needed, initial_rels);
		if (enable_geqo && levels_needed >= geqo_threshold)
			return geqo(root, levels_needed, initial_rels);
		return standard_join_search(root, levels_needed, initial_rels);
	}

	result = dphyp(root, levels_needed, initial_rels);

	/* Successfully found join order */
	if (list_length(result) == 1)
	{
		RelOptInfo *rel = linitial(result);
		list_free(result);
		return rel;
	}

	/* 
	 * Single relation in List means we successfully found query plan.
	 * But we may fail and in this case we can have:
	 * 
	 * 1. NIL
	 * 2. Any amount of relations
	 * 
	 * First case means we can not do anything, so pass 'initial_rels' to conventional DPsize/GEQO.
	 * The second case means we have implicit joins, but plans for disjoint sets are found - pass what we have found to DPsize/GEQO.
	 */
	if (result == NIL)
	{
		result = initial_rels;
	}

	if (prev_join_search_hook)
		return prev_join_search_hook(root, list_length(result), result);
	if (enable_geqo && list_length(result) >= geqo_threshold)
		return geqo(root, list_length(result), result);
	return standard_join_search(root, list_length(result), result);
}

void
_PG_init(void)
{

	DefineCustomBoolVariable("pg_dphyp.enabled",
							 "pg_dphyp join enumeration algorithm is enabled",
							 NULL,
							 &enabled,
							 enabled,
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
