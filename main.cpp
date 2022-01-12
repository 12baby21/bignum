#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
#include <limits.h>

    using namespace std;

extern int P;

/**
 * invert module
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
**/
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

int main()
{
    bn op1(11);
    bn op2(31);
    int len1 = 32;
    int len2 = 32;

    int len = 1;
    while (len < len1 * 2 || len < len2 * 2)
        len <<= 1;
    ntt(&op1, len, 1);
    ntt(&op2, len, 1);
    bn ans;
    for (int i = 0; i < len; i++)
    {
        ans.array[i] = op1.array[i] * op2.array[i] % P;
    }
    /**
    这里需要考虑溢出的问题
    for (int i = 0; i < len; i++)
    {
        ans.array[i] += op1.array[i];
        if (ans.array[i] > INT_MAX)
        {
            ans.array[i + 1] += ans.array[i] / (INT_MAX+1);
            ans.array[i] %= 10;
        }
    }

    int pos = 0;
    for (int i = len - 1; i >= 0; i--)
    {
        if (ans[i])
        {
            pos = i;
            break;
        }
    }
    **/

    ntt(&ans, 32, true);
    bn_print(&ans);

    return 0;
}
