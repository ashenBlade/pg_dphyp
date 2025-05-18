#ifndef SIMPLEBMS_H
#define SIMPLEBMS_H

#include "nodes/bitmapset.h"

bitmapword bmw_make_singleton(int x);
bitmapword bmw_all_bit_set(int x);
bitmapword bmw_add_member(bitmapword a, int x);
bitmapword bmw_next_member(bitmapword a, int prevbit);
bitmapword bmw_prev_member(bitmapword a, int prevbit);
bitmapword bmw_union(bitmapword a, bitmapword b);
bitmapword bmw_difference(bitmapword a, bitmapword b);
bitmapword bmw_is_member(bitmapword a, int x);
bitmapword bmw_overlap(bitmapword a, bitmapword b);
bool bmw_is_subset(bitmapword a, bitmapword b);
bool bmw_equal(bitmapword a, bitmapword b);
int bmw_first(bitmapword a);

static inline bool bmw_is_empty(bitmapword bmw)
{
	return bmw == 0;
}

bool bmw_single_element(bitmapword bmw, int x);

uint32 bmw_hash_value(bitmapword bmw);
uint32 bmw_hash(const void *key, Size keysize);
int bmw_match(const void *key1, const void *key2, Size keysize);

#endif
