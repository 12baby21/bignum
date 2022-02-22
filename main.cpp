#include <iostream>
#include "bn.h"
#include "MontG.h"
#include "time.h"
#include "common.h"
#include <limits.h>
#include <cmath>
#include <algorithm>

using namespace std;


int main()
{
    // This function is used to test s1ome atomic operations
    /* 1. Fermat inverse */
/*
    bn p(13);
    for (int i = 7; i <= 89; i++)
    {
        bn a(i);
        bn inverse;
        bn_Fermat_inverse(&inverse, &a, &p);

        bn_print(&inverse);
    }
*/

    /* 2. montMult in MontG.cpp */
/*
    bn x(10);
    bn y(7);
    bn m(51);
    // bn r2m(256);
    bn res1;
    montMult(&x, &y, &m, 6, &res1);
    bn res2;
    bn_mul(&res2, &x, &y);
    bn_mod(&res2, &res2, &m);
    // modExp(&x, &e, 3, &m, 4, &r2m, &res);
    bn_print(&res1);
    cout << "gold:" << endl;
    bn_print(&res2);
*/

    /* 3. rshift and lshift */
/*    
    bn x;
    x.array[0] = 1;
    cout << "Origin number:" << endl;
    bn_print(&x);

    bn_rshift_bits(&x, &x, 700);
    cout << "number after rshift" << endl;
    bn_print(&x);

    bn_lshift_bits(&x, &x, 10);
    cout << "number after lshift" << endl;
    bn_print(&x);
*/    

    /* 4. Take low k bits */
/*
    bn x;
    x.array[0] = 6554;
    x.array[1] = 257485;
    x.array[2] = 58430;
    x.array[3] = 354689;
    bn_print(&x);
    bn_TakeLowBits(&x, &x, 104);
    bn_print(&x);
*/

    /* 5. barrett reduction */
    bn res;
    bn a(146119744);
    bn b(19789);
    BarrettReduction(&res, &a, &b);
    bn_print(&res);

    /* 6. bn_qmod */
/*
    bn res;
    bn a(51666);
    bn b(2);
    bn c(19789);
    bn_qmod(&res, &a, &b, &c);
    bn_print(&res);
*/
    return 0;
}


/*
int main()
{
    bn m(1);
    bn c;
    bn res;

    bn n, g, lambda, mu, n_2;
    GenKey(1024, &n, &g, &lambda, &mu, &n_2);
    cout<<"n:"<<endl;
    bn_print(&n);
    cout<<endl;
    cout<<"g:"<<endl;
    bn_print(&g);
    cout<<endl;
    cout<<"lambda:"<<endl;
    bn_print(&lambda);
    cout<<endl;
    cout<<"mu:"<<endl;
    bn_print(&mu);
    cout<<endl;
    cout<<"n^2:"<<endl;
    bn_print(&n_2);
    cout<<endl;



    Encryption(&c, &m, &g, &n, &n_2);
    Decryption(&res, &c, &lambda, &n, &n_2);

    cout<<"m:"<<endl;
    bn_print(&m);
    cout<<endl;

    cout<<"res:"<<endl;
    bn_print(&res);

    return 0;
}
*/

/**
int main()
{
    clock_t start = clock();

    bn RandNum(551541242518964165);
    //ProduceRandom(&RandNum);
    RandNum.array[0] |= 1;

    cout << "We have generated a random number:" << endl;
    bn_print(&RandNum);

    bn p;
    bn_nextprime(&p, &RandNum);
    cout << "The next prime of RandNum is: ";
    bn_print(&p);

    clock_t finish = clock();
    cout << "Time elapsed:" << dec << 1.0 * (finish - start) / CLOCKS_PER_SEC << "s" << endl;
    return 0;
}
**/
