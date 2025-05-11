#ifndef DPHYPSIMPLE_H
#define DPHYPSIMPLE_H

#include "nodes/pathnodes.h"

List *dphyp_simple(PlannerInfo *root, int levels_needed, List *initial_rels);

#endif