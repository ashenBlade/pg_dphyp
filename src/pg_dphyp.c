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

#include "unionset.h"
#include "dphyp_generic.h"
#include "dphyp_simple.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

static join_search_hook_type prev_join_search_hook = NULL;

static bool enabled = true;

void _PG_init(void);
void _PG_fini(void);

static RelOptInfo *dphyp_join_search(PlannerInfo *root, int levels_needed, List *initial_rels)
{
	List *result;

	if (!enabled)
	{
		if (prev_join_search_hook)
			return prev_join_search_hook(root, levels_needed, initial_rels);
		if (enable_geqo && levels_needed >= geqo_threshold)
			return geqo(root, levels_needed, initial_rels);
		return standard_join_search(root, levels_needed, initial_rels);
	}

	if (BITS_PER_BITMAPWORD < levels_needed)
		result = dphyp_generic(root, levels_needed, initial_rels);
	else
		result = dphyp_simple(root, levels_needed, initial_rels);

	Assert(result != NIL);
	if (list_length(result) == 1)
		return linitial(result);

	elog(NOTICE, "running implicit join path");
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
