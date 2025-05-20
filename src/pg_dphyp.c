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

#include "dphyp.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

static join_search_hook_type prev_join_search_hook = NULL;

/* GUC */
/* Extension is enabled and should run DPhyp */
static bool dphyp_enabled = true;
/* Do not apply DPhyp if SJ found (LEFT/RIGHT/OUTER etc...) */
static bool dphyp_skip_sj = true;
/* Threshold for GEQO used by DPhyp */
static int dphyp_geqo_threshold = 12;

void _PG_init(void);
void _PG_fini(void);

static inline bool contains_sj(PlannerInfo *root)
{
	/* 
	 * There may be more fine-grained checks, but for small OLTP queries
	 * this will introduce too much overhead.
	 * So, just check whether or not we have SJ in query at all.
	 */
	return root->join_info_list != NIL;
}

static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	List *result;
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

	result = dphyp(root, levels_needed, initial_rels);	

	
	/* Successfully found join order */
	if (list_length(result) == 1)
	{
		RelOptInfo *rel = linitial(result);
		list_free(saved_join_rel_list);
		list_free(result);
		return rel;
	}
	
	/* Restore state before proceeding to DPsize/GEQO */
	list_free(root->join_rel_list);
	root->join_rel_list = saved_join_rel_list;

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
		if (prev_join_search_hook)
			return prev_join_search_hook(root, levels_needed, initial_rels);
		if (enable_geqo && levels_needed >= geqo_threshold)
			return geqo(root, levels_needed, initial_rels);
		return standard_join_search(root, levels_needed, initial_rels);
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
