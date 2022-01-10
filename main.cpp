#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
using namespace std;

/**
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
int main()
{
    bn RandNum(14292);
    ProduceRandom(&RandNum);
    clock_t start = clock();

    cout << "We have generated a random number:" << endl;
    bn_print(&RandNum);

    bn p;
    bn_nextprime(&p, &RandNum);
    //cout << rabinmiller(&RandNum, 25)<<endl;
    cout << "The next prime of RandNum is: ";
    bn_print(&p);

    clock_t finish = clock();
    cout << "Time elapsed:" << dec << finish - start << "us" << endl;
    return 0;
}