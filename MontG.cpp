#include <iostream>
#include <time.h>
#include "bn.h"
#include "MontG.h"

using namespace std;

#define LUT_SIZE 2049

static bn lut[LUT_SIZE];

/*
 * Performs bitwise Montgomery modular multiplication ( X*Y*R^(-k) mod M),k=# of bits in x or y
 *
 * Parameters:
 * 		x,y,m - bignums
 * 		mBits - # of bits in m
 * 		out	  - bignum result
 */

void montMult(bn_ptr x, bn_ptr y, bn_ptr m, int mBits, bn_ptr out)
{
	bignum t;
	int i;

	int y0Bit = y->getSize(); //-------y0bit仅需要算一次
	
	for(i = 0; i < mBits ; i++)
    {// efficient loop exit

		// 计算公式 t' <-- (Y * X_i + t + op * m) * R^(-1)
		int t0Bit = bn_getbit(&t, 0);
		int xiBit = bn_getbit(x, i);	//loop exit requires subtraction here,
		
		// op <-- (Y_0 * X_i + t_0)* S mod R ; 
        // (S = -M^(-1) mod R ， S = 1 ,   R为进制数 R = 2)
		int op = t0Bit + (xiBit & y0Bit);  

		if(xiBit == 1)
        {// t' <-- (Y * X_i + t)
			bn_add(&t, y, &t);
		}
	
		if(op == 1)
        {// t' <-- (Y * X_i + t + op * m)
			bn_add(&t, m, &t);
		}

		//t' <-- (Y * X_i + t + q * n) * R^(-1)
		_bn_rshift(&t,&t, 1);
	}
	// 此时 t < 2M
	if(bn_cmp(&t, m) == LARGER){
		bn_sub(&t, m, &t);
	}

	bn_assign(out, &t);
}



/*mod exp, no LUT */
void modExp(bn_ptr x, 
            bn_ptr e, 
            int eBits, 
            bn_ptr m, 
            int mBits, 
            bn_ptr r2m,
            bn_ptr out)
{

	bn z, one;
	bn parr[3];
	bn zarr[3];

	// ρ = r^(2*K), 即输入参数 r2m
	// zarr ← MontMul(1,ρ) ,进入蒙哥马利域：out = 1 * r^k mod n
	bn_assign(&z, 1);
	montMult(&z,r2m,m, mBits, &zarr[1]);

	// parr ← MontMul(x,ρ), 进入蒙哥马利域：out = x * r^k mod n
	montMult(x,r2m,m, mBits,&parr[1]);

	bignum tm;
	// 从e的低位开始扫描计算， L-R二进制算法 (left to right Binary Algorithm)?
	// 由于平方与相乘互相独立，因此可以并行计算，
	// 但由于e的每1bit不可能都为1，因此倍的资源不能达到两倍的性能提升
	int i = 0;
	for(; i < eBits; i++){

		bn_assign(&tm, &parr[1]);
		montMult(&tm,&parr[1],m, mBits, &parr[2]);  // 平方

		if(bn_getbit(e, i) == 1){ // 相乘
			montMult(&zarr[1],&parr[1],m,mBits,&zarr[2]);
		}else{
			bn_assign(&zarr[2],&zarr[1]);
		}

		//printf("num bits p: %d, num bits z: %d\n", bignum_numbits(&parr[1]), bignum_numbits(&zarr[1]));
		// 这里的数组不需要移位，只需要加条件判断让数组交替使用即可
		bn_assign(&parr[1], &parr[2]);
		bn_assign(&zarr[1], &zarr[2]);
	}

	bn_assign(&one, 1);
	montMult(&zarr[1], &one, m, mBits, out);
}


/* Mod Exp using precomputed LUT */
void modExpLUT(bn_ptr x, bn_ptr e, int eBits, bn_ptr m, int mBits, bn_ptr r2m, bn_ptr out)
{
	bignum z,one;
	bignum parr[3];
	bignum zarr[3];

	//reduce z?
	bn_assign(&z, 1);
	montMult(&z, r2m, m, mBits, &zarr[1]);

	//MontMul(x,ρ) 可以提前计算好
	bn_assign(&parr[1], &lut[0]);

	// int b = 1;
	int i = 0;
	for(; i < eBits; i++){
		bn_assign(&parr[2],&lut[i+1]);

		if(bn_getbit(e, i) == 1){
			montMult(&zarr[1], &lut[i], m, mBits, &zarr[2]);
		}else{
			bn_assign(&zarr[2], &zarr[1]);
		}

		//printf("num bits p: %d, num bits z: %d\n", bignum_numbits(&parr[1]), bignum_numbits(&zarr[1]));
		bn_assign(&parr[1], &parr[2]);
		bn_assign(&zarr[1], &zarr[2]);
		// b++;
	}

	bn_assign(&one, 1);
	montMult(&zarr[1], &one, m, mBits, out);
}

int parseLUT(int start, int num){

	FILE *f;
	f = fopen("LUT.txt", "r");
	char str[num+1];
	static char* zpad[8] = {"", "0", "00", "000", "0000", "00000", "000000", "0000000"};
	//printf("Parse LUT: %d to %d\n", start, num + start);

	//puts("");
	int i = 0;
	int a = 0;
	for(; i < (num + start); i++){
		fscanf(f, "%s\n", str);

		if(i >= start){
			int len = strlen(str);
			char* dup;
			if((len & 1) == 1){
				sprintf(str, "%s%s", zpad[1], (dup = strdup(str)));
				len++;
			}

			int lenMod8 = len -((len >> 3) << 3);
			if(lenMod8 != 0){
				sprintf(str, "%s%s", zpad[lenMod8], (dup = strdup(str)));
				len += lenMod8;
			}
			bignum_from_string(&lut[i-start], str, len);
		}
		a++;
	}
	return i;
}

