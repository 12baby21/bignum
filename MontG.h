#ifndef __MONTG_H__
#define __MONTG_H__
#include "bn.h"


void modExp(bn_ptr x, bn_ptr e, int eBits, bn_ptr m, int mBits, bn_ptr r2m,  bn_ptr out);
void montMult(bn_ptr x, bn_ptr y, bn_ptr m, int mBits, bn_ptr out);


#endif