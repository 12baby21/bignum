#include <iostream>
#include "bn.h"
#include "time.h"
using namespace std;

int main()
{
    bignum a(1024);
    bignum b(5);
    bignum c(57);
    bignum res;

    bignum d(1024);
    bignum e(5);
    bignum f(57);
    bignum res2;

    clock_t start_time = clock();
    bn_qmod(&res, &a, &b, &c);
    clock_t finish_time = clock();
    bn_print(&res);
    cout << "Time elapsed: " << dec << \
        finish_time - start_time << "ms" << endl;

    start_time = clock();
    bn_pow(&res2, &d, &e);
    bn_mod(&res2, &res2, &f);
    finish_time = clock();
    bn_print(&res2);
    cout << "Time elapsed: " << dec << \
        finish_time - start_time << "ms" << endl;
    return 0;
}