/** TODO List
 **	NTT multiplication
 ** paillier
 *
 **/
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include "bn.h"
#include "common.h"
#include "primes.h"
using namespace std;
/* Functions for shifting number in-place. */
static void _lshift_one_bit(bignum *a);
static void _rshift_one_bit(bignum *a);
static void _lshift_word(bignum *a, int nwords);
static void _rshift_word(bignum *a, int nwords);

bignum::bignum()
{
	int i = 0;
	for (i = 0; i < BN_ARRAY_SIZE; i++)
	{
		array[i] = 0;
	}
}

bignum::bignum(uint64_t num)
{
	// bignum()
	int i = 0;
	for (i = 0; i < BN_ARRAY_SIZE; i++)
	{
		array[i] = 0;
	}
#ifdef WORD_SIZE
#if (WORD_SIZE == 1)
	// low 8 bits -> high 8 bits
	array[0] = (num & 0x000000ff);
	array[1] = (num & 0x0000ff00) >> 8;
	array[2] = (num & 0x00ff0000) >> 16;
	array[3] = (num & 0xff000000) >> 24;
#elif (WORD_SIZE == 2)
	// low 16 bits -> high 16 bits
	array[0] = (num & 0x0000ffff);
	array[1] = (num & 0xffff0000) >> 16;
#elif (WORD_SIZE == 4)
	array[0] = num;
	DTYPE_TMP num_32 = 32;
	DTYPE_TMP tmp = num >> num_32; /* bit-shift with U64 operands to force 64-bit results */
	array[1] = tmp;
#endif
#endif
}

int bignum::getSize()
{
	int size = BN_ARRAY_SIZE;
	for (int i = BN_ARRAY_SIZE - 1; i >= 0; --i)
	{
		if (this->array[i] != 0)
			return size;
	}
	return 0;
}

/* often used bignum */
bn one(1);
bn zero(0);
bn two(2);
bn three(3);
bn five(5);
bn seven(7);

/* ntt used Primitive root */
const int M = 2013265921;
const int g = 31;
LL mul_x1, mul_x2, mul_res;
const int UnitMax = 32768;

/**	will be used if data struct is specified
#if(WORD_SIZE == 1)
#define UnitMax = 256
#elif(WORD_SIZE == 2)
#define UnitMax = 32768
#elif(WORD_SIZE == 4)
#define UnitMax = 65536
#endif
**/
void bn_assign(bignum *op1, bignum *op2)
{
	for (int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		op1->array[i] = op2->array[i];
	}
}

void bn_assign(bignum *op1, DTYPE_TMP n)
{
	/* Endianness issue if machine is not little-endian? */
#ifdef WORD_SIZE
#if (WORD_SIZE == 1)
	op1->array[0] = (n & 0x000000ff);
	op1->array[1] = (n & 0x0000ff00) >> 8;
	op1->array[2] = (n & 0x00ff0000) >> 16;
	op1->array[3] = (n & 0xff000000) >> 24;
#elif (WORD_SIZE == 2)
	op1->array[0] = (n & 0x0000ffff);
	op1->array[1] = (n & 0xffff0000) >> 16;
#elif (WORD_SIZE == 4)
	op1->array[0] = n & MAX_VAL;
	op1->array[1] = n >> 32;
#endif
#endif
}

/* Basic arithmetic operations for bignum: */
void bn_add(bignum *res, bignum *op1, bignum *op2)
{
	DTYPE_TMP tmp;
	int carry = 0;
	for (int i = 0; i < BN_ARRAY_SIZE; ++i)
	{
		tmp = (DTYPE_TMP)op1->array[i] + op2->array[i] + carry;
		carry = (tmp > MAX_VAL);
		res->array[i] = (tmp & MAX_VAL);
	}
}

void bn_sub(bignum *res, bignum *op1, bignum *op2)
{
	DTYPE_TMP tmp1;
	DTYPE_TMP tmp2;
	DTYPE_TMP inter_res;
	bn tmp_res;
	int borrow = 0;
	for (int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		tmp1 = (DTYPE_TMP)op1->array[i] + (MAX_VAL + 1);
		tmp2 = (DTYPE_TMP)op2->array[i] + borrow;
		inter_res = (tmp1 - tmp2);
		tmp_res.array[i] = (DTYPE)(inter_res & MAX_VAL);
		borrow = (inter_res <= MAX_VAL);
	}
	bn_assign(res, &tmp_res);
}

void bn_mul(bignum *res, bignum *op1, bignum *op2)
{
	int i, j;

	DTYPE lowbits = 0;
	DTYPE highbits = 0; // carry
	bignum tmp_res;
	for (i = 0; i < BN_ARRAY_SIZE; ++i)
	{
		bignum row;
		for (j = 0; j < BN_ARRAY_SIZE; ++j)
		{
			if (i + j < BN_ARRAY_SIZE)
			{
				bignum tmp;
				DTYPE_TMP intermediate = ((DTYPE_TMP)op1->array[i] * (DTYPE_TMP)op2->array[j]);
				lowbits = (intermediate & MAX_VAL);
				highbits = (intermediate >> 32);

				tmp.array[0] = lowbits;
				tmp.array[1] = highbits;

				_lshift_word(&tmp, i + j);
				bn_add(&row, &tmp, &row);
			}
		}
		bn_add(&tmp_res, &row, &tmp_res);
	}
	bn_assign(res, &tmp_res);
}

void bn_div(bignum *res, bignum *dividend, bignum *divisor)
{
	bignum current(1);
	bignum denom; // fenmu
	bignum tmp;
	bignum tmp_res;
	bn_assign(&denom, divisor);
	bn_assign(&tmp, dividend);

	const DTYPE_TMP half_max = 1 + (DTYPE_TMP)(MAX_VAL / 2);
	bool overflow = false;
	while (bn_cmp(&denom, dividend) != LARGER) // while(denom <= dividend)
	{
		if (denom.array[BN_ARRAY_SIZE - 1] >= half_max)
		{
			overflow = true;
			break;
		}
		_lshift_one_bit(&current);
		_lshift_one_bit(&denom);
	}
	if (!overflow)
	{
		_rshift_one_bit(&denom);   // denom >>= 1;
		_rshift_one_bit(&current); // current >>= 1;
	}

	while (!bn_is_zero(&current)) // while (current != 0)
	{
		if (bn_cmp(&tmp, &denom) != SMALLER) //   if (dividend >= denom)
		{
			bn_sub(&tmp, &tmp, &denom);			 //     dividend -= denom;
			bn_or(&tmp_res, &current, &tmp_res); //     answer |= current;
		}
		_rshift_one_bit(&current); //   current >> 1;
		_rshift_one_bit(&denom);   //   denom >> 1;
	}							   // return answer;
	bn_assign(res, &tmp_res);
}

void bn_pow(bignum *res, bignum *base, bignum *power)
{
	if (bn_cmp(power, &zero) == EQUAL)
	{
		/* Return 1 when exponent is 0 -- n^0 = 1 */
		bn_inc(res);
	}
	else
	{
		bignum tmp;
		bn_assign(res, base);
		bn_dec(power);

		/* Begin summing products: */
		while (!bn_is_zero(power))
		{
			/* res = res * base */
			bn_mul(res, res, base);
			/* Decrement power by one */
			bn_dec(power);
		}
	}
}

void bn_mod(bignum *res, bignum *a, bignum *b)
{
	bignum mulcb;
	bignum temp_res;
	bignum c;

	/* c = (a / b) */
	bn_div(&c, a, b);

	/* mulcb = (c * b) */
	bn_mul(&mulcb, &c, b);

	/* c = a - tmp */
	bn_sub(&temp_res, a, &mulcb);

	bn_assign(res, &temp_res);
}

void bn_qmod(bignum *res, bignum *a, bignum *b, bignum *c)
{
	bn tmp_res(1);
	for (int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		DTYPE pow = b->array[i];
		if (pow == 0)
		{
			for (int j = 0; j < WORD_SIZE * 8; j++)
			{
				bn_mul(a, a, a);
				bn_mod(a, a, c);
			}
		}
		else
		{
			while (pow)
			{
				if (pow & 1)
				{
					bn_mul(&tmp_res, &tmp_res, a);
					bn_mod(&tmp_res, &tmp_res, c);
				}
				bn_mul(a, a, a);
				bn_mod(a, a, c);
				pow = pow >> 1;
			}
		}
	}
	bn_assign(res, &tmp_res);
}

void bn_qpow(bn_ptr res, bn_ptr base, bn_ptr power)
{
	bn tmp_res(1);
	for (int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		DTYPE pow = power->array[i];
		if (pow == 0)
		{
			for (int j = 0; j < WORD_SIZE * 8; j++)
			{
				bn_mul(base, base, base);
			}
		}
		else
		{
			while (pow)
			{
				if (pow & 1)
				{
					bn_mul(&tmp_res, &tmp_res, base);
				}
				bn_mul(base, base, base);
				pow = pow >> 1;
			}
		}
	}
	bn_assign(res, &tmp_res);
}

void bn_inc(bignum *n)
{
	bn_add(n, n, &one);
}

void bn_dec(bignum *n)
{
	bn_sub(n, n, &one);
}


// 费马小定理求逆元
// inverse = a ^ {p-2}
void bn_Fermat_inverse(bn_ptr inverse, bn_ptr a, bn_ptr p)
{
	// p-2
	bn p_2;
	bn_sub(&p_2, p, &two);

	bn_qpow(inverse, a, &p_2);
}


/* optimized operator */
static void ButterflyOperation(int* x, int len)
{
	int i, j, k;
	for (i = 1, j = len >> 1; i < len - 1; i++)
	{
		if (i < j)
			swap(x[i], x[j]);
		k = len >> 1;
		while (j >= k)
		{
			j -= k;
			k /= 2;
		}
		if (j < k)
			j += k;
	}
}

static void NTT(int* x, int len, bool isINTT = false)
{
	/**
	 * 涉及的运算
	 * 模幂 -> qmod
	 * 模乘 -> * %
	 */
	ButterflyOperation(x, len);
	for (int h = 2; h <= len; h <<= 1)
	{
		//
		int gn = qmod(g, (M - 1) / h, M);
		if (isINTT)
		{
			gn = qmod(gn, M - 2, M);
		}
		for (int j = 0; j < len; j += h)
		{
			int g0 = 1;
			// 这是一个计算域
			for (int k = j; k < j + h / 2; k++)
			{
				LL u = x[k];
				LL t = (LL)g0 * x[k + h / 2] % M;
				x[k] = (u + t) % M;
				x[k + h / 2] = (u - t + M) % M;
				g0 = (LL)g0 * gn % M;
			}
		}
	}
	if (isINTT == true)
	{
		int inv = qmod(len, M - 2, M);

		for (int i = 0; i < len; i++)
		{
			LL tmp = x[i];
			tmp = tmp * inv % M;
			x[i] = tmp;
		}
	}
}

static void preProcessOperand(bn_ptr op, int* buffer)
{
	for(int i = 0; i < op->actual_array_size; ++i)
	{
		*(buffer + i) = op->array[i];
	}
}

// 乘法域
// 对于输入的乘数，需要预先读取操作数至更大的有限域中
void bn_ntt_mul(bn_ptr res, int &n, bn_ptr op1, int n1, bn_ptr op2, int n2)
{
	
	n = 1;
	int _n = n1 + n2 - 1;
	// fill "len" to 2^l
	while (n < _n)
		n <<= 1;
	// coefficient -> point value

	// 把操作数存在int型有限域中
	// 也许可以当作队列
	int buffer1[BN_ARRAY_SIZE];
	int buffer2[BN_ARRAY_SIZE];
	int res_buffer[BN_ARRAY_SIZE];

	// 可以并行处理
	// 把操作数预先读入乘法域中的buffer
	preProcessOperand(op1, buffer1);
	preProcessOperand(op2, buffer2);
	
	// 系数表达式 -> 点值表达式
	NTT(buffer1, n, false);
	NTT(buffer2, n, false);


	// dot multiplication
	for (int i = 0; i < n; ++i)
	{
		mul_x1 = buffer1[i];
		mul_x2 = buffer2[i];
		mul_res = mul_x1 * mul_x2 % M;
		res_buffer[i] = mul_res;
	}
	// 系数表达式 -> 点值表达式
	NTT(res_buffer, n, true);

	// 解决乘法溢出	SSA
	for (int i = 0; i < n; i++)
	{
		// 在硬件电路中直接取高位为进位，低位为保留位即可
		res->array[i + 1] +=res_buffer[i] / UnitMax;
		res->array[i] = res_buffer[i] % UnitMax;
	}
	while (res->array[n - 1] == 0)
	{
		n--;
	}
	res->actual_array_size = n;
}

/* Basic arithmetic operations for basic types: */
int qmod(int a, int b, int M)
{
	int ans = 1;
	while (b)
	{
		if (b & 1)
			ans = (LL)ans * a % M;
		a = (LL)a * a % M;
		b >>= 1;
	}
	return ans;
}

/* Bitwise operations: */
void bn_or(bignum *res, bignum *op1, bignum *op2)
{
	int i;
	for (i = 0; i < BN_ARRAY_SIZE; ++i)
	{
		res->array[i] = (op1->array[i] | op2->array[i]);
	}
}

void bn_and(bignum *res, bignum *op1, bignum *op2)
{

	int i;
	for (i = 0; i < BN_ARRAY_SIZE; ++i)
	{
		res->array[i] = (op1->array[i] & op2->array[i]);
	}
}

/* Special operators and comparison */
int bn_cmp(bignum *op1, bignum *op2)
{
	int i = BN_ARRAY_SIZE;
	do
	{
		i -= 1; // do decrement to start with the last element
		if (op1->array[i] > op2->array[i])
		{
			return 1;
		}
		else if (op1->array[i] < op2->array[i])
		{
			return -1;
		}
	} while (i != 0);
	return 0;
}

void bn_invert(bn_ptr lambdainvert, bn_ptr lambda, bn_ptr n)
{
	bn gcds, gcdt, gcdgcd;
	extended_euclidean(lambdainvert, &gcdt, &gcdgcd, lambda, n);
}

int bn_getbit(const bn_ptr a, int n)
{
	int index = n / (WORD_SIZE * 8);
	int dst = n - index * (WORD_SIZE * 8);
	return (a->array[index] >> dst) & 1;
}

void bn_setbit(const bn_ptr a, int n)
{
	int array_index = n / 32;
	int sub_index = n % 32;
	int i = 1;
	while (n--)
		i = i << 1;
	a->array[array_index] |= i;
}

void _bn_rshift(bn_ptr a, bn_ptr b, int nbits)
{
	/* Handle shift in multiples of word-size */
	int nwords = nbits >> 5;
	if (nwords != 0)
	{
		int z = nwords << 5;
		_rshift_word(a, nwords);
		nbits -= (z);
	}

	if (nbits != 0)
	{
		int z = 32 - nbits;
		int i;
		for (i = 0; i < (BN_ARRAY_SIZE - 1); i++)
		{
			// a->array[i + 1] << (z) 即取高位的低32-nbits位
			DTYPE higher = a->array[i+1] << z;
			DTYPE lower = a->array[i] >> nbits;
			a->array[i] = higher | lower;
		}
		a->array[i] >>= nbits;
	}
}


int bn_numbits(bignum *bn)
{
	int n = BN_ARRAY_SIZE - 1;
	int b;
	for (; n >= 0; n--)
	{
		b = bn_getbit(bn, n);
		if (b == 1)
		{
			return n + 1;
		}
	}
	return 0;
}

bool bn_is_zero(bignum *n)
{
	int i;
	for (i = 0; i < BN_ARRAY_SIZE; i++)
	{
		if (n->array[i])
		{
			return 0;
		}
	}
	return 1;
}

void bn_print(bignum *num)
{
	for (int i = BN_ARRAY_SIZE - 1; i >= 0; i--)
	{
		cout << hex << num->array[i] << " ";
	}
	cout << endl;
}

/* Functions for shifting number in-place. */
static void _lshift_one_bit(bignum *a)
{
	require(a, "a is null");

	int i;
	for (i = (BN_ARRAY_SIZE - 1); i > 0; --i)
	{
		a->array[i] = (a->array[i] << 1) | (a->array[i - 1] >> ((8 * WORD_SIZE) - 1));
	}
	a->array[0] <<= 1;
}

static void _rshift_one_bit(bignum *a)
{
	require(a, "a is null");

	int i;
	for (i = 0; i < (BN_ARRAY_SIZE - 1); ++i)
	{
		a->array[i] = (a->array[i] >> 1) | (a->array[i + 1] << ((8 * WORD_SIZE) - 1));
	}
	a->array[BN_ARRAY_SIZE - 1] >>= 1;
}

static void _lshift_word(bignum *a, int nwords)
{
	int i;
	/* Shift whole words */
	for (i = (BN_ARRAY_SIZE - 1); i >= nwords; --i)
	{
		a->array[i] = a->array[i - nwords];
	}
	/* Zero pad shifted words. */
	for (; i >= 0; --i)
	{
		a->array[i] = 0;
	}
}

static void _rshift_word(bignum *a, int nwords)
{
	register int i = i - i;
	for (; i < nwords; i += 4)
	{
		a->array[i] = a->array[i + 1];
		a->array[i + 1] = a->array[i + 2];
		a->array[i + 2] = a->array[i + 3];
		a->array[i + 3] = a->array[i + 4];
	}
	register int z = z - z;
	for (; i < BN_ARRAY_SIZE; i += 4)
	{
		a->array[i] = z;
		a->array[i + 1] = z;
		a->array[i + 2] = z;
		a->array[i + 3] = z;
	}
}

/* functions for big primes */
int calc_trial_divisions(int bits)
{
	if (bits <= 512)
		return 64;
	else if (bits <= 1024)
		return 128;
	else if (bits <= 2048)
		return 384;
	else if (bits <= 4096)
		return 1024;
	return 2048;
}

void ProduceRandom(bn_ptr RandNum)
{
	clock_t t;
	unsigned int RandomNumber; //记录随机数

	for (int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		t = clock();
		srand((unsigned)t); // seed
		RandomNumber = rand();
		RandNum->array[i] = RandomNumber;
	}
	// mask
	RandNum->array[BN_ARRAY_SIZE - 1] |= 0x80000000;
}

bool rabinmiller(bn_ptr n, int trails)
{
	int s = 0;
	bn temp;
	bn n_1; // n-1
	bn_sub(&temp, n, &one);
	bn_sub(&n_1, n, &one);

	// let n-1 present in the way of (2^s)*t
	int thisBit = temp.array[0] & 0x1;
	int isZero = bn_is_zero(&temp);
	// isZero can be accelerated
	while (thisBit == 0 && !isZero)
	{
		_rshift_one_bit(&temp);
		s++;
		thisBit = temp.array[0] & 0x1;
		isZero = bn_is_zero(&temp);
	}
	bn t;
	bn_assign(&t, &temp);

	// trials times Miller-Rabin test
	while (trails--)
	{
		bignum b(primes[100 - trails]);
		bignum y;

		bn_qmod(&y, &b, &t, n);

		int cmp1 = bn_cmp(&y, &one);
		int cmp2 = bn_cmp(&y, &n_1);
		if (cmp1 == 0 || cmp2 == 0)
			return true;
		// y != n-1
		for (int j = 1; j <= (s - 1) && bn_cmp(&y, &n_1) != 0; ++j)
		{
			bn_qmod(&y, &y, &two, n);
			if (bn_cmp(&y, &one) == 0)
				return false;
		}
		if (bn_cmp(&y, &n_1) != 0)
			return false;
	}
	return true;
}

void bn_nextprime(bn_ptr p, bn_ptr n)
{
	unsigned long difference;
	int i;
	unsigned prime_limit;
	unsigned long prime;
	unsigned incr;

	/* First handle tiny numbers */
	// if n < 2, p = 2 is the next prime
	if (bn_cmp(n, &two) < 0)
	{
		bn_assign(p, &two);
		return;
	}

	// if n >= 2, set p as an odd larger or equal to n
	bn_add(p, n, &one);
	bn_setbit(p, 0);

	// if p <= 7, i.e. p = 5 or 7, then p is a prime
	if (bn_cmp(p, &seven) <= 0)
		return;

	prime_limit = 166;

	/* Compute residues modulo small odd primes */
	bn moduli[166];
	for (;;)
	{
		// prime[1]
		for (i = 0; i < prime_limit; i++)
		{
			bn prime(primes[i + 1]);
			bn_mod(moduli + i, p, &prime);
		}

#define INCR_LIMIT 0x10000 /* deep science */
		bn diff;
		for (difference = incr = 0; incr < INCR_LIMIT; difference += 2)
		{
			for (i = 0; i < prime_limit; i++)
			{
				bn prime(primes[i + 1]);
				bn r;
				bn inc(incr);
				bn sum;
				bn_add(&sum, moduli + i, &inc);
				bn_mod(&r, &sum, &prime);

				// r == 0 means that (moduli + inc) is a composite number
				// i.e. (p + inc) is a composite number
				if (bn_is_zero(&r))
					goto next;
			}

			{
				bn tmp(difference);
				bn_assign(&diff, &tmp);
				bn_add(p, p, &diff);
				difference = 0;

				/* Miller-Rabin test */
				// if test passed, then p is probable a prime
				if (rabinmiller(p, 25))
					goto done;
			}
		next:;
			incr += 2;
		}
		bn_add(p, p, &diff);
		bn_print(p);
		difference = 0;
	}
done:;
}
