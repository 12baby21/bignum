#ifndef COMMON_H
#define COMMON_H
#include "bn.h"



void GenKey(uint64_t bits, bn_ptr n, bn_ptr g, bn_ptr lambda, bn_ptr mu, bn_ptr n_2);
void Encryption(bn_ptr c, bn_ptr m, bn_ptr g, bn_ptr n, bn_ptr n_2);
void Decryption(bn_ptr res, bn_ptr c, bn_ptr lambda, bn_ptr n, bn_ptr n_2);



#endif