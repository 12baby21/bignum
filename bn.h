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

#ifndef BITS
  #define BITS 1024
#endif


#define BN_ARRAY_SIZE ((BITS / WORD_SIZE / 8))


/* Here comes the compile-time specialization for how large the underlying array size should be. */
/* The choices are 1, 2 and 4 bytes in size with uint32, uint64 for WORD_SIZE==4, as temporary. */
#ifndef WORD_SIZE
  #error Must define WORD_SIZE to be 1, 2, 4
#elif (WORD_SIZE == 1)
  /* Data type of array in structure */
  #define DTYPE                    uint8_t
  /* bitmask for getting MSB */
  #define DTYPE_MSB                ((DTYPE_TMP)(0x80))
  /* Data-type larger than DTYPE, for holding intermediate results of calculations */
  #define DTYPE_TMP                uint32_t
  /* sprintf format string */
  #define SPRINTF_FORMAT_STR       "%.02x"
  #define SSCANF_FORMAT_STR        "%2hhx"
  /* Max value of integer type */
  #define MAX_VAL                  ((DTYPE_TMP)0xFF)
#elif (WORD_SIZE == 2)
  #define DTYPE                    uint16_t
  #define DTYPE_TMP                uint32_t
  #define DTYPE_MSB                ((DTYPE_TMP)(0x8000))
  #define SPRINTF_FORMAT_STR       "%.04x"
  #define SSCANF_FORMAT_STR        "%4hx"
  #define MAX_VAL                  ((DTYPE_TMP)0xFFFF)
#elif (WORD_SIZE == 4)
  #define DTYPE                    uint32_t
  #define DTYPE_TMP                uint64_t
  #define DTYPE_MSB                ((DTYPE_TMP)(0x80000000))
  #define SPRINTF_FORMAT_STR       "%.08x"
  #define SSCANF_FORMAT_STR        "%8x"
  #define MAX_VAL                  ((DTYPE_TMP)0xFFFFFFFF)
#endif
#ifndef DTYPE
  #error DTYPE must be defined to uint8_t, uint16_t uint32_t or whatever
#endif

/* Custom assert macro - easy to disable */
#define require(p, msg) assert(p && #msg)




class bignum
{
public:
	DTYPE array[BN_ARRAY_SIZE];

	bignum();
	bignum(uint64_t num);
};

typedef bignum* bn_ptr;
typedef bignum bn;


/* Tokens returned by bignum_cmp() for value comparison */
enum { SMALLER = -1, EQUAL = 0, LARGER = 1 };

/* Initialization functions: */
/* op1 = op2 */
void bn_assign(bignum *op1, bignum *op2);   
/* op1 = n */
void bn_assign(bignum *op1, DTYPE_TMP n);       

/* Basic arithmetic operations: */
/* res = op1 + op2 */
void bn_add(bignum* res, bignum* op1, bignum* op2); 
/* res = op1 - op2 */
void bn_sub(bignum* res, bignum* op1, bignum* op2); 
/* res = op1 * op2 */
void bn_mul(bignum* res, bignum* op1, bignum* op2); 
/* res = dividend / divisor */
void bn_div(bignum* res, bignum* dividend, bignum* divisor); 
/* res = base ^ power */
void bn_pow(bignum* res, bignum* base, bignum* power);       
/* res = a mod b */
void bn_mod(bignum* res, bignum* a, bignum* b);
/* res = a^b mod c O(logN) but hasn't theoretical acceleration */
void bn_qmod(bignum* res, bignum* a, bignum* b, bignum* c);


/* Bitwise operations: */
/* res = op1 | op2 */
void bn_or(bignum* res, bignum* op1, bignum* op2);  
/* res = op1 & op2 */
void bn_and(bignum* res, bignum* op1, bignum* op2);

/* Special operators and comparison */
/* Compare: returns LARGER=1, EQUAL=0 or SMALLER=-1 */
int bn_cmp(bignum *op1, bignum *op2);   
/* For comparison with zero */            
bool bn_is_zero(bignum* n);   
/* Increment: add one to n */                      
void bn_inc(bignum* n);    
/* Decrement: sub one from n */                        
void bn_dec(bignum* n); 
/* invert */
void bn_invert(bn_ptr lambdainvert, bn_ptr lambda, bn_ptr n);


/* Utility Functions */
/* get the n-th bit of a */
int bn_getbit(const bignum* a, int n);
/* set the n-th bit to 1 */
void bn_setbit(const bn_ptr a, int n);
/* return the bits number of bn */
int bn_numbits(bignum *bn);               
/* print the bignum in hex type */
void bn_print(bignum* num);               


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


