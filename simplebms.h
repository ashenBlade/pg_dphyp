/* 
 * Helpful set of functions that operate on single bitmapword
 * instead of whole Bitmapset that must be allocated on heap.
 * DPhyp heavily depends on efficient bit operations, so operating
 * on single bitmapword seems good performance optimization.
 */

#ifndef SIMPLEBMS_H
#define SIMPLEBMS_H

#include "nodes/bitmapset.h"
#include "port/pg_bitutils.h"

#if PG_MAJORVERSION_NUM < 17
#if BITS_PER_BITMAPWORD == 32
#define bmw_leftmost_one_pos(w)			pg_leftmost_one_pos32(w)
#define bmw_rightmost_one_pos(w)		pg_rightmost_one_pos32(w)
#else
#define bmw_leftmost_one_pos(w)			pg_leftmost_one_pos64(w)
#define bmw_rightmost_one_pos(w)		pg_rightmost_one_pos64(w)
#endif
#endif

#if PG_MAJORVERSION_NUM < 15
#define pg_rotate_left32(x,k)	(((x)<<(k)) | ((x)>>(32-(k))))
#endif

#define ValidateBmwPosition(x) Assert(0 <= (x) && (x) < BITS_PER_BITMAPWORD)
#define MAKE_BMW(x)  ((bitmapword) (1 << (x)))

/* 
 * Create bitmapword with single bit set at 'x' position
 */
static inline bitmapword
bmw_make_singleton(int x)
{
	ValidateBmwPosition(x);
	return MAKE_BMW(x);
}

/* 
 * Create bitmapword with all bits prior to 'x' position set
 */
static inline bitmapword
bmw_all_bit_set(int x)
{
	ValidateBmwPosition(x);
	return ((bitmapword)(1 << (x + 1))) - 1;
}

/* 
 * Add member 'x' to 'bmw' bitmapword
 */
static inline bitmapword
bmw_add_member(bitmapword bmw, int x)
{
	ValidateBmwPosition(x);
	return bmw | MAKE_BMW(x);
}

/* 
 * Get all elements from 'a' without elements from 'b'
 */
static inline bitmapword
bmw_difference(bitmapword a, bitmapword b)
{
	return a & ~b;
}

/* 
 * Check that 'a' is subset of 'b'
 */
static inline bool
bmw_is_subset(bitmapword a, bitmapword b)
{
	return (a & b) == a;
}

/* 
 * Check if 'x' is member of 'bmw'
 */
static inline bitmapword
bmw_is_member(bitmapword bmw, int x)
{
	ValidateBmwPosition(x);
	return (bmw & MAKE_BMW(x)) != 0;
}

/* 
 * Check if 'a' and 'b' have any common members
 */
static inline bool
bmw_overlap(bitmapword a, bitmapword b)
{
	return (a & b) != 0;
}

/* 
 * Get index of first member of bitmapword from start.
 * Used to get representative of hypernode.
 */
static inline int
bmw_first(bitmapword bmw)
{
	if (bmw == 0)
		return 0;
	else
		return bmw_rightmost_one_pos(bmw);
}

static inline int
bmw_lowest_bit(bitmapword bmw)
{
	return bmw & (-bmw);
}

/* 
 * Get next member of 'bmw' starting from 'prevbit'.
 * Pass -1 to 'prevbit' at the start.
 * Returns -1 if there are no more members.
 */
static inline bitmapword
bmw_next_member(bitmapword bmw, int prevbit)
{
	bitmapword mask;
	
	if (prevbit != -1)
		ValidateBmwPosition(prevbit);

	mask = (~(bitmapword) 0) << (prevbit + 1);
	bmw &= mask;

	if (bmw == 0)
		return -1;
	
	return bmw_rightmost_one_pos(bmw);
}

/* 
 * Get previous member of 'bmw' starting from 'prevbit'.
 * Pass -1 to 'prevbit' at the start.
 * Returns -1 if there are no more members.
 */
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
		ValidateBmwPosition(prevbit);
		--prevbit;
	}

	mask = (~(bitmapword) 0) >> (BITS_PER_BITMAPWORD - (prevbit + 1));
	bmw &= mask;

	if (bmw == 0)
		return -1;
	
	return bmw_leftmost_one_pos(bmw);
}

/* 
 * Hash function for bitmapword to be used in HTAB
 */
static inline bitmapword
bmw_hash_value(bitmapword x)
{
	/* Copied from hashfn.c */
	uint32		a,
				b,
				c;

	a = b = c = 0x9e3779b9 + (uint32) sizeof(uint32) + 3923095;
	a += (int32) x;

	c ^= b; c -= pg_rotate_left32(b,14);
	a ^= c; a -= pg_rotate_left32(c,11);
	b ^= a; b -= pg_rotate_left32(a,25);
	c ^= b; c -= pg_rotate_left32(b,16);
	a ^= c; a -= pg_rotate_left32(c, 4);
	b ^= a; b -= pg_rotate_left32(a,14);
	c ^= b; c -= pg_rotate_left32(b,24);

#if BITS_PER_BITMAPWORD == 64

	a += (int32) (x >> 32);

	c ^= b; c -= pg_rotate_left32(b,14);
	a ^= c; a -= pg_rotate_left32(c,11);
	b ^= a; b -= pg_rotate_left32(a,25);
	c ^= b; c -= pg_rotate_left32(b,16);
	a ^= c; a -= pg_rotate_left32(c, 4);
	b ^= a; b -= pg_rotate_left32(a,14);
	c ^= b; c -= pg_rotate_left32(b,24);

#endif

	/* report the result */
	return c;
}

/* 
 * Generic hash function for bitmapword
 */
static inline uint32
bmw_hash(const void *key, Size keysize)
{
	Assert(keysize == sizeof(bitmapword));

	return bmw_hash_value(*(const bitmapword *)key);
}

/* 
 * Comparison function for bitmapword members in HTAB
 */
static inline int
bmw_match(const void *key1, const void *key2, Size keysize)
{
	Assert(keysize == sizeof(bitmapword));
	
	return *(const bitmapword *)key1 != *(const bitmapword *)key2;
}

/* 
 * Check that 'bmw' contains only single member 'x'
 */
static inline bool
bmw_single_element(bitmapword bmw, int x)
{
	ValidateBmwPosition(x);
	return bmw == MAKE_BMW(x);
}

/* 
 * Check that 'bmw' has only single bit set.
 * Does not check that 'bmw' is empty.
 */
static inline bool
bmw_is_singleton(bitmapword bmw)
{
	return (bmw & (bmw - 1)) == 0;
}


/* 
 * Check if 'bmw' is empty
 */
static inline bool
bmw_is_empty(bitmapword bmw)
{
	return bmw == 0;
}

#endif
