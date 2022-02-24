#ifndef BN_H
#define BN_H
#include <stdio.h>

/*

Big number library - arithmetic on multiple-precision unsigned integers.

This library is an implementation of arithmetic on arbitrarily large integers.

The difference between this and other implementations, is that the data structure
has optimal memory utilization (i.e. a 1024 bit integer takes up 128 bytes RAM),
and all memory is allocated statically: no dynamic allocation for better or worse.

Primary goals are correctness, clarity of code and clean, portable implementation.
Secondary goal is a memory footprint small enough to make it suitable for use in
embedded applications.


The current state is correct functionality and adequate performance.
There may well be room for performance-optimizations and improvements.

*/

#include <stdint.h>
#include <assert.h>

/* This macro defines the word size in bytes of the array that constitutes the big-number data structure. */
// default: 1 byte 32 bits
#ifndef WORD_SIZE
#define WORD_SIZE 4
#endif

// length of key
#ifndef BITS
#define BITS 1024
#endif

// 数组的每个单元存储WORD_SIZE位数据
// 2倍的存储空间为乘法服务
#define BN_ARRAY_SIZE ((BITS / WORD_SIZE / 8) * 2)

/* Here comes the compile-time specialization for how large the underlying array size should be. */
/* The choices are 1, 2 and 4 bytes in size with uint32, uint64 for WORD_SIZE==4, as temporary. */
#ifndef WORD_SIZE
#error Must define WORD_SIZE to be 1, 2, 4
// 1 byte
#elif (WORD_SIZE == 1)
/* Data type of array in structure */
#define DTYPE uint8_t
#define DTYPE_MSB ((DTYPE_TMP)(0x80))
#define DTYPE_TMP uint32_t
#define MAX_VAL ((DTYPE_TMP)0xFF)
// 2 bytes
#elif (WORD_SIZE == 2)
#define DTYPE uint16_t
#define DTYPE_TMP uint32_t
#define DTYPE_MSB ((DTYPE_TMP)(0x8000))
#define MAX_VAL ((DTYPE_TMP)0xFFFF)
// int
#elif (WORD_SIZE == 4)
#define DTYPE uint32_t
#define DTYPE_TMP uint64_t  // used for multiplication
#define DTYPE_MSB ((DTYPE_TMP)(0x80000000))   // used for fetch the maximum bit
#define MAX_VAL ((DTYPE_TMP)0xFFFFFFFF)
#endif
#ifndef DTYPE
#error DTYPE must be defined to uint8_t, uint16_t uint32_t or whatever
#endif

class bignum
{
public:
  DTYPE array[BN_ARRAY_SIZE];
  bignum();
  bignum(uint64_t num);
  bignum(bignum* num);
  int32_t size;   // size < 0 represents negative, size > 0 represents positive
  int32_t bitnum; // bit number of 
  // get the actual size of the bignum
  int getSize();
  int getBitnum();
};

typedef bignum *bn_ptr;
typedef bignum bn;
typedef long long LL;

/* Tokens returned by bignum_cmp() for value comparison */
enum
{
  SMALLER = -1,
  EQUAL = 0,
  LARGER = 1
};

/* Initialization functions: */
/* op1 = op2 */
void bn_assign(bignum *op1, bignum *op2);
/* op1 = n */
void bn_assign(bignum *op1, DTYPE_TMP n);

/* Basic arithmetic operations for bignum: */
/* res = op1 + op2 */
void bn_add(bignum *res, bignum *op1, bignum *op2);
/* res = op1 - op2 */
void bn_sub(bignum *res, bignum *op1, bignum *op2);
/* res = op1 * op2 */
void bn_mul(bignum *res, bignum *op1, bignum *op2);
/* res = dividend / divisor */
void bn_div(bignum *res, bignum *dividend, bignum *divisor);
/* res = base ^ power */
void bn_pow(bignum *res, bignum *base, bignum *power);
/* res = a mod b */
void bn_mod(bignum *res, bignum *a, bignum *b);
/* res = a^b mod c O(logN) but hasn't theoretical acceleration */
void bn_qmod(bignum *res, bignum *a, bignum *b, bignum *c);
/* res = base ^ power */
void bn_qpow(bn_ptr res, bn_ptr base, bn_ptr power);
/* res = x ^ n, this function is for cpu-calculable operands */
LL quick(LL x, LL n);

/* Basic arithmetic operations for basic types: */
/* return: a^b mod M */
int qmod(int a, int b, int M);

/* optimized operator */
void bn_ntt_mul(bn_ptr res, int &n, bn_ptr op1, int n1, bn_ptr op2, int n2);
/* res = a mod b, note that bit number of a must less than 2 * bit number of b */
void BarrettReduction(bn_ptr res, bn_ptr a, bn_ptr b);


/* Bitwise operations: */
/* res = op1 | op2 */
void bn_or(bignum *res, bignum *op1, bignum *op2);
/* res = op1 & op2 */
void bn_and(bignum *res, bignum *op1, bignum *op2);
/* a = a >> n, inplace */
void _bn_rshift(bn_ptr a, bn_ptr b, int nbits);
/* res = a >> nbits */
void bn_rshift_bits(bn_ptr res, bn_ptr a, int nbits);
/* res = a << nbits */
void bn_lshift_bits(bn_ptr res, bn_ptr a, int nbits);
/* res = a << nwords */
void bn_lshift_words(bn_ptr res, bn_ptr a, int nwords);
/* res = a >> nwords */
void bn_rshift_words(bn_ptr res, bn_ptr a, int nwords);
/* res =  a[0, highbit)*/
void bn_TakeLowBits(bn_ptr res, bn_ptr a, int highbit);

/* Special operators and comparison */
/* Compare: returns LARGER=1, EQUAL=0 or SMALLER=-1 */
int bn_cmp(bignum *op1, bignum *op2);
/* For comparison with zero */
bool bn_is_zero(bignum *n);
/* Increment: add one to n */
void bn_inc(bignum *n);
/* Decrement: sub one from n */
void bn_dec(bignum *n);
/* Fermat inverse */
void bn_Fermat_inverse(bn_ptr inverse, bn_ptr a, bn_ptr p);

/* Utility Functions */
/* get the n-th bit of a */
int bn_getbit(const bn_ptr a, int n);
/* set the n-th bit to 1 */
void bn_setbit(const bn_ptr a, int n);
/* return the bits number of bn */
int bn_numbits(bignum *bn);
/* print the bignum in hex type */
void bn_print(bn_ptr num);
/* get a range of bits */

/* functions for big primes */
/* Calculate trail number */
int calc_trial_divisions(int bits);
/* Generate a random odd number */
void ProduceRandom(bn_ptr RandNum);
/* Miller-Rabin Prime Test */
bool rabinmiller(bn_ptr n, int trails);
/* set p as the next prime of n */
void bn_nextprime(bn_ptr p, bn_ptr n);

#endif /* #ifndef __BIGNUM_H__ */
