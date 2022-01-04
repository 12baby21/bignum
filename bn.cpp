#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "bn.h"
#include "common.h"
using namespace std;
/* Functions for shifting number in-place. */
static void _lshift_one_bit(bignum* a);
static void _rshift_one_bit(bignum* a);
static void _lshift_word(bignum* a, int nwords);
static void _rshift_word(bignum* a, int nwords);


bignum::bignum()
{
	int i = 0;
	for(i = 0; i < BN_ARRAY_SIZE; i++)
	{
		array[i] = 0;
	}
}

bignum::bignum(uint64_t num)
{
	// bignum()
	int i = 0;
	for(i = 0; i < BN_ARRAY_SIZE; i++)
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

bignum one(1);
bignum zero(0);
bignum two(2);



void bn_assign(bignum *op1, bignum *op2)
{
	for(int i = 0; i < BN_ARRAY_SIZE; i++)
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


/* Basic arithmetic operations: */
void bn_add(bignum* res, bignum* op1, bignum* op2)
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

void bn_sub(bignum* res, bignum* op1, bignum* op2)
{
	DTYPE_TMP tmp1;
	DTYPE_TMP tmp2;
	DTYPE_TMP inter_res;
	int borrow = 0;
	for(int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		tmp1 = (DTYPE_TMP)op1->array[i] + (MAX_VAL + 1);
		tmp2 = (DTYPE_TMP)op2->array[i] + borrow;
		inter_res = (tmp1 - tmp2);
		res->array[i] = (DTYPE)(inter_res & MAX_VAL);
		borrow = (inter_res <= MAX_VAL);
	}
}

void bn_mul(bignum* res, bignum* op1, bignum* op2)
{
    int i, j;

	DTYPE lowbits = 0;
	DTYPE highbits = 0;		// carry
	bignum tmp_res;			// avoid if op1|op2 = res
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

void bn_div(bignum* res, bignum* dividend, bignum* divisor)
{
	bignum current(1);
	bignum denom;		// fenmu
	bignum tmp;
	bn_assign(&denom, divisor);
	bn_assign(&tmp, dividend);

	const DTYPE_TMP half_max = 1 + (DTYPE_TMP)(MAX_VAL / 2);
	bool overflow = false;
	while(bn_cmp(&denom, dividend) != LARGER)	// while(denom <= dividend)
	{
		if(denom.array[BN_ARRAY_SIZE -1] >= half_max)
		{
			overflow = true;
			break;
		}
		_lshift_one_bit(&current);
		_lshift_one_bit(&denom);
	}
	if(!overflow)
	{
		_rshift_one_bit(&denom);                  // denom >>= 1;
		_rshift_one_bit(&current);                // current >>= 1;
	}


	while (!bn_is_zero(&current))           // while (current != 0)
	{
		if (bn_cmp(&tmp, &denom) != SMALLER)  //   if (dividend >= denom)
		{
		bn_sub(&tmp, &tmp, &denom);         //     dividend -= denom;
		bn_or(res, &current, res);              //     answer |= current;
		}
		_rshift_one_bit(&current);                //   current >> 1;
		_rshift_one_bit(&denom);                  //   denom >> 1;
	}                                           // return answer;
}

void bn_pow(bignum* res, bignum* base, bignum* power)
{
	if(bn_cmp(power, &zero) == EQUAL)
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

void bn_mod(bignum* res, bignum* a, bignum* b)
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

void bn_qmod(bignum* res, bignum* a, bignum* b, bignum* c)
{
	bn_assign(res, 1);
	for(int i = 0; i < BN_ARRAY_SIZE; i++)
	{
		DTYPE pow = b->array[i];
		while(pow)
		{
			if(pow & 1)
			{
				bn_mul(res, res, a);
				bn_mod(res, res, c);
			}				
			bn_mul(a, a, a);
			bn_mod(a, a, c);
			pow = pow >> 1;
		}
	}
}

void bn_inc(bignum* n)
{
  bn_add(n, n, &one);
}

void bn_dec(bignum* n)
{
	bn_sub(n, n, &one);
}


/* Bitwise operations: */
void bn_or(bignum* res, bignum* op1, bignum* op2)
{
  int i;
  for (i = 0; i < BN_ARRAY_SIZE; ++i)
  {
    res->array[i] = (op1->array[i] | op2->array[i]);
  }
}

void bn_and(bignum* res, bignum* op1, bignum* op2)
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
		i -= 1;	// do decrement to start with the last element
		if(op1->array[i] > op2->array[i])
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


int bn_getbit(const bignum* a, int n)
{
	int index = n / (WORD_SIZE * 8);
	int dst = n - index * (WORD_SIZE * 8);
	return (a->array[index] >> dst) & 1;
}

int bn_numbits(bignum *bn)
{
	int n = BN_ARRAY_SIZE -1;
	int b;
	for (; n >= 0; n--){
		b = bn_getbit(bn, n);
		if (b == 1){
			return n+1;
		}
	}
	return 0;
}

bool bn_is_zero(bignum* n)
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

void bn_print(bignum* num){
    for(int i = BN_ARRAY_SIZE-1; i >= 0; i--)
    {
        cout << hex << num->array[i] << " ";
    }
	cout << endl;
}


/* Functions for shifting number in-place. */
static void _lshift_one_bit(bignum* a)
{
  require(a, "a is null");

  int i;
  for (i = (BN_ARRAY_SIZE - 1); i > 0; --i)
  {
    a->array[i] = (a->array[i] << 1) | (a->array[i - 1] >> ((8 * WORD_SIZE) - 1));
  }
  a->array[0] <<= 1;
}

static void _rshift_one_bit(bignum* a)
{
  require(a, "a is null");

  int i;
  for (i = 0; i < (BN_ARRAY_SIZE - 1); ++i)
  {
    a->array[i] = (a->array[i] >> 1) | (a->array[i + 1] << ((8 * WORD_SIZE) - 1));
  }
  a->array[BN_ARRAY_SIZE - 1] >>= 1;
}

static void _lshift_word(bignum* a, int nwords)
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

static void _rshift_word(bignum* a, int nwords)
{
	register int i = i-i;
	for (; i < nwords; i+=4)
	{
		a->array[i]   = a->array[i + 1];
		a->array[i+1] = a->array[i + 2];
		a->array[i+2] = a->array[i + 3];
		a->array[i+3] = a->array[i + 4];
	}
	register int z = z-z;
	for (; i < BN_ARRAY_SIZE; i+=4)
	{
		a->array[i]   = z;
		a->array[i+1] = z;
		a->array[i+2] = z;
		a->array[i+3] = z;
	}
}


/* functions for big primes */

/*使用三个rand()生成伪随机数组合生成一个奇数随机数，作为伪素数 
**系统时间为种子 
**并返回生成的这个大奇数 
*/ 
unsigned int ProduceRandomOdd(){
    //UINT无符号整形，各伪随机数放在RandomArray数组中 
    time_t t;//c++时间类型 
    unsigned int RandomNumber;//记录随机数 
	do{
    
		srand((unsigned)time(&t));//srand(seed)用于给rand()函数设定种子,此处用系统时间
		//生成 
        RandomNumber=(rand()<<17)|(rand()<<3)|(rand()); 
        //cout<<RandomNumber<<endl;
	}while(RandomNumber%2==0||RandomNumber<100000000); 
	//返回   
	return RandomNumber;
}

//模重复平方算法求(b^n)%m
size_t repeatMod(size_t base, size_t n, size_t mod){
    
    size_t a = 1;
    while(n){
    
        if(n&1){
    
            a=(a*base)%mod;
        }
        base=(base*base)%mod;
        n=n>>1;
    }
    return a;
}

//Miller-Rabin素数检测
bool rabinmiller(size_t n, size_t k){
    
	
    int s=0;
    int temp=n-1;    
	
	//将n-1表示为(2^s)*t  
    while ((temp&0x1)==0&&temp){
    
        temp=temp>>1;
        s++;
    }   
    size_t t = temp;
    
	//判断k轮误判概率不大于(1/4)^k
    while(k--){
    
        srand((unsigned)time(0));
        size_t b = rand()%(n-2)+2; //生成一个b(2≤a ≤n-2)

        size_t y = repeatMod(b,t,n); 
        if (y == 1 || y == (n-1))
            return true;
        for(int j = 1; j<=(s-1) && y != (n-1); ++j){
    
            y = repeatMod(y,2,n);
            if (y == 1)
                return false;
        }
        if (y != (n-1))
            return false;
    }
    return true;
}

/**

bool is_probable_prime(bignum *num, int trials)
{
	int isLarger = bn_cmp(num, &two);
	assert(isLarger == 1);
	if(isLarger == 0)
		return true;
	bignum residual;
	bn_mod(&residual, num, &two);

	int s = 0;
	bignum d;
	bn_sub(&d, num, &one);
	// To present n - 1 as 2^s * d
	while(trails)


}

**/