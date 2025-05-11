#include "postgres.h"
#include "utils/palloc.h"
#include "nodes/bitmapset.h"

#include "unionset.h"

void us_init(unionset_state *state, int size)
{
	int *leaders;
	leaders = palloc(sizeof(int) * size);
	for (size_t i = 0; i < size; i++)
	{
		leaders[i] = i;
	}

	state->leaders = leaders;
	state->size = size;
}

int us_leader(unionset_state *state, int node)
{
	Assert(node < state->size);
	if (state->leaders[node] == node)
		return node;
	else
		return state->leaders[node] = us_leader(state, state->leaders[node]);
}

void us_union(unionset_state *state, int a, int b)
{
	int a_leader;
	int b_leader;

	a_leader = us_leader(state, a);
	b_leader = us_leader(state, b);

	state->leaders[a_leader] = b_leader;
}

List *us_collect(unionset_state *state)
{
	List *result;
	Bitmapset **disjoint_sets;
	
	disjoint_sets = palloc0(sizeof(Bitmapset *) * state->size);
	for (size_t i = 0; i < state->size; i++)
	{
		int leader = us_leader(state, i);
		if (disjoint_sets[leader] == NULL)
			disjoint_sets[leader] = bms_make_singleton(i);
		else
			disjoint_sets[leader] = bms_add_member(disjoint_sets[leader], i);
	}

	result = NIL;
	for (size_t i = 0; i < state->size; i++)
	{
		if (disjoint_sets[i] == NULL)
			continue;
		
		result = lappend(result, disjoint_sets[i]);
	}
	
	pfree(disjoint_sets);
	return result;
}

void us_free(unionset_state *state)
{
	pfree(state->leaders);
}
