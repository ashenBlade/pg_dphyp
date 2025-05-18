#ifndef SIMPLEBMS_H
#define SIMPLEBMS_H

#include "nodes/bitmapset.h"
#include "port/pg_bitutils.h"
#include "common/hashfn.h"

#define ValidateBmwMember(x) Assert(0 <= (x) && (x) < BITS_PER_BITMAPWORD)
#define MAKE_BMW(x)  ((bitmapword) (1 << (x)))

static inline bitmapword
bmw_make_singleton(int x)
{
	ValidateBmwMember(x);
	return MAKE_BMW(x);
}

static inline bitmapword
bmw_all_bit_set(int x)
{
	ValidateBmwMember(x);
	return ~((bitmapword)0) >> (BITS_PER_BITMAPWORD - (x + 1));
}

static inline bitmapword
bmw_add_member(bitmapword bmw, int x)
{
	ValidateBmwMember(x);
	return bmw | MAKE_BMW(x);
}

static inline bitmapword
bmw_union(bitmapword a, bitmapword b)
{
	return a | b;
}

/* 
 * All elements from 'a' without elements from 'b'
 */
static inline bitmapword
bmw_difference(bitmapword a, bitmapword b)
{
	return a & ~b;
}

/* 
 * 'a' is subset of 'b'
 */
static inline bool
bmw_is_subset(bitmapword a, bitmapword b)
{
	return (a | b) == b;
}

static inline bool
bmw_equal(bitmapword a, bitmapword b)
{
	return a == b;
}

static inline bitmapword
bmw_is_member(bitmapword bmw, int x)
{
	ValidateBmwMember(x);
	return (bmw & MAKE_BMW(x)) != 0;
}

static inline bitmapword
bmw_overlap(bitmapword a, bitmapword b)
{
	return (a & b) != 0;
}

static inline int
bmw_first(bitmapword bmw)
{
	if (bmw == 0)
		return -1;
	else
		return bmw_rightmost_one_pos(bmw);
}

static inline bitmapword
bmw_next_member(bitmapword bmw, int prevbit)
{
	bitmapword mask;
	
	if (prevbit != -1)
		ValidateBmwMember(prevbit);

	mask = (~(bitmapword) 0) << (prevbit + 1);
	bmw &= mask;

	if (bmw == 0)
		return -1;
	
	return bmw_rightmost_one_pos(bmw);
}

static inline bitmapword
bmw_prev_member(bitmapword bmw, int prevbit)
{
	bitmapword mask;

	if (prevbit == 0)
		return -1;

	if (prevbit == -1)
	{
		prevbit = BITS_PER_BITMAPWORD - 1;
	}
	else
	{
		ValidateBmwMember(prevbit);
		--prevbit;
	}

	mask = (~(bitmapword) 0) >> (BITS_PER_BITMAPWORD - (prevbit + 1));
	bmw &= mask;

	if (bmw == 0)
		return -1;
	
	return bmw_leftmost_one_pos(bmw);
}

static inline uint32
bmw_hash_value(bitmapword bmw)
{
#if BITS_PER_BITMAPWORD == 32
	return bmw;
#else /* BITS_PER_BITMAPWORD == 64 */
	return hash_bytes_uint32((uint32)bmw) ^ hash_bytes_uint32((uint32)(bmw >> 32));
#endif
}

static inline uint32
bmw_hash(const void *key, Size keysize)
{
	Assert(keysize == sizeof(bitmapword));
	return bmw_hash_value(*(const bitmapword *)key);
}

static inline int
bmw_match(const void *key1, const void *key2, Size keysize)
{
	Assert(keysize == sizeof(bitmapword));
	
	return *(const bitmapword *)key1 != *(const bitmapword *)key2;
}

static inline bool
bmw_single_element(bitmapword bmw, int x)
{
	return pg_popcount((const char *)&bmw, sizeof(bmw)) == 1 && bmw_is_member(bmw, x);
}

static inline bool
bmw_is_empty(bitmapword bmw)
{
	return bmw == 0;
}

#endif
