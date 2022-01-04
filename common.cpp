#include <iostream>
#include "bn.h"
#include "common.h"
using namespace std;



void extended_euclidean(bn_ptr gcds, bn_ptr gcdt, bn_ptr gcdgcd, bn_ptr a, bn_ptr b)
{
	if (b == 0) 
    { //b=0時直接結束求解
		bn_assign(gcds, (uint64_t)1);
        bn_assign(gcdt, (uint64_t)0);
        bn_assign(gcdgcd, (uint64_t)0);
        return;
	}
	
    bn old_r;
    bn_assign(&old_r, a);
    bn r;
    bn_assign(&r, b);
	bn old_s(1);
    bn s(0);
	bn old_t(0);
    bn t(1);
    bn zero(0);
	while (bn_cmp(&r, &zero) != EQUAL) 
    { //按擴展歐基里德算法進行循環
		bn q;
        bn_div(&q, &old_r, &r);
		bn temp;
        bn_assign(&temp, &old_r);
		bn_assign(&old_r, &r);
        bn mulqr;
        bn_mul(&mulqr, &q, &r);
		bn_sub(&r, &temp, &mulqr);
		bn_assign(&temp , &old_s);
		bn_assign(&old_s, &s);
        bn mulqs;
        bn_mul(&mulqs, &q, &s);
        bn_assign(&s, (uint64_t)0);
		bn_sub(&s, &temp, &mulqs);
		bn_assign(&temp, &old_t);
		bn_assign(&old_t, &t);
        bn mulqt;
        bn_mul(&mulqt, &q, &t);
        bn_sub(&t, &temp, &mulqt);
	}
    bn_assign(gcds, &old_s);
    bn_assign(gcdt, &old_t);
    bn_assign(gcdgcd, &old_r);
}


