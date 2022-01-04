#include <iostream>
#include "bn.h"

using namespace std;

bignum yi(1);

void GenKey(uint64_t bits, bignum* n, bignum* g, bignum* lambda, bignum* mu, bignum* n_2)
{

    bignum PrimeP(17);
    bignum PrimeQ(19);

    // Generate Public Key (n, g)
    // n = PrimeP * PrimeQ
    // g = n + 1, if PrimeP and PrimeQ have the same length
    bn_mul(n, &PrimeP, &PrimeQ);
    bn_add(g, n, &yi);

    // PrimeP = PrimeP-1
    // PrimeQ = PrimeQ-1
    bn_sub(&PrimeP, &PrimeP, &yi);
    bn_sub(&PrimeQ, &PrimeQ, &yi);

    // Generate Secret Key (lambda, mu)
    // lambda = (PrimeP-1) * (PrimeQ-1)
    // mu = lambda -1
    bn_mul(lambda, &PrimeP, &PrimeQ);
    bn_sub(mu, lambda, &yi);

    // n_2 = n^2
    bn_mul(n_2, n, n);
}

void Encryption(bn_ptr c, bn_ptr m, bn_ptr g, bn_ptr n, bn_ptr n_2)
{
    // Random Key
    bignum r(2);

    // gm = g^m mod n^2
    // rn = r^n mod n^2
    bignum gm;
    bignum rn;

    // c = g^m * r^n mod n^2 = (g^m mod n^2) * (r^n mod n^2) mod n^2
    // c = gm * rn mod n^2
    bn_qmod(&gm, g, m, n_2);
    bn_qmod(&rn, &r, n, n_2);
    bn_mul(c, &gm, &rn);
    bn_mod(c, c, n_2);
}

void Decryption(bn_ptr res, bn_ptr c, bn_ptr lambda, bn_ptr n, bn_ptr n_2)
{
    bn l;
    cout<<"c = g^m * r^n mod n^2:"<<endl;
    bn_print(c);
    bn_qmod(&l, c, lambda, n_2);
    bn_sub(&l, &l, &yi);
    bn_div(&l, &l, n);

    bn lambdainvert;
    bn_invert(&lambdainvert, lambda, n);
    bn_mul(&l, &l, &lambdainvert);
    bn_mod(res, &l, n);

}