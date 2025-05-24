#ifndef DPHYPSIMPLE_H
#define DPHYPSIMPLE_H

#include "nodes/pathnodes.h"

/*
 * Entry point for DPhyp join ordering algorithm.
 * Returns a list that should be interpreted basing on it's size:
 * 
 * 0  - DPhyp can not find valid query plan, so it should pass control to standard planner
 * 1  - query plan found and single element is required plan
 * >1 - there are implicit joins found, so each element is a plan for such disjoint graph and this list should be passed to standard planner
 */
RelOptInfo *dphyp(PlannerInfo *root, int levels_needed, List *initial_rels);

#endif