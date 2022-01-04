#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
using namespace std;


int main()
{
    bn a(1);
    bn c;
    bn m;
    bn res;

    bn n, g, lambda, mu, n_2;
    GenKey(1024, &n, &g, &lambda, &mu, &n_2);
    Encryption(&c, &m, &g, &n, &n_2);
    Decryption(&res, &c, &lambda, &n, &n_2);

    bn_print(&a);
    cout<<endl;
    bn_print(&c);
    cout<<endl;
    bn_print(&res);

    return 0;
}

