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
    bn RandNum;
    ProduceRandomOdd(&RandNum);
    bool flag = false;
    clock_t start = clock();
    while (!flag)
    {
        cout << "We have generated a random number:" << endl;
        bn_print(&RandNum);
        int trails = calc_trial_divisions(1024);

        bool isPrime = rabinmiller(&RandNum, trails);

        cout << "This number ";
        if (isPrime)
        {
            flag = true;
            cout << "is ";
        }
        else
        {
            bn two(2);
            cout << "isn't ";
            bn_add(&RandNum, &RandNum, &two);
        }
        cout << "a prime number!" << endl;
    }
    clock_t finish = clock();
    cout << "Time elapsed:" << dec << finish - start << "us" << endl;
    return 0;
}