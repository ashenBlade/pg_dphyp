#ifndef UNIONSET_H
#define UNIONSET_H

#include "nodes/pg_list.h"

typedef struct unionset_state
{
	int *leaders;
	int size;
} unionset_state;

void us_init(unionset_state *state, int size);
int us_leader(unionset_state *state, int node);
void us_union(unionset_state *state, int a, int b);
List *us_collect(unionset_state *state);
void us_free(unionset_state *state);

#endif