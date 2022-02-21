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
    bn r(2);

    // gm = g^m mod n^2
    // rn = r^n mod n^2
    bn gm;
    bn rn;
    // c = g^m * r^n mod n^2 = (g^m mod n^2) * (r^n mod n^2) mod n^2
    // c = gm * rn mod n^2
    bn_qmod(&gm, g, m, n_2);
    bn_qmod(&rn, &r, n, n_2);
    bn_mul(c, &gm, &rn);
    bn_mod(c, c, n_2);
    cout << "g^m:" << endl;
    bn_print(&gm);
    cout << "r^n:" << endl;
    bn_print(&rn);
}

void Decryption(bn_ptr res, bn_ptr c, bn_ptr lambda, bn_ptr n, bn_ptr n_2)
{
    bn l;
    cout<<"c = g^m * r^n mod n^2:"<<endl;
    bn_print(c);
    cout<<endl;
    cout<<"c^lambda mod n^2:"<<endl;
    bn_qmod(&l, c, lambda, n_2);
    bn_print(&l);
    cout<<endl;
    bn_sub(&l, &l, &yi);
    cout<<"l = c^lambda mod n^2-1:"<<endl;
    bn_print(&l);
    cout<<endl;
    bignum l1_n;    // (l-1) / n
    bn_div(&l1_n, &l, n);
    cout<<"l = (l-1)/n:"<<endl;
    bn_print(&l1_n);
    cout<<endl;

    bn lambdainvert;
    bn_Fermat_inverse(&lambdainvert, lambda, n);
    bn_mul(&l, &l, &lambdainvert);
    bn_mod(res, &l, n);

}