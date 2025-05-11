#ifndef DPHYP_GENERIC_H
#define DPHYP_GENERIC_H

#include "nodes/pathnodes.h"

/* 
 * Возвращает список непересекающихся RelOptInfo *
 */
List *dphyp_generic(PlannerInfo *root, int levels_needed, List *initial_rels);

#endif