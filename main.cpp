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
    bn RandNum(327);
    //ProduceRandomOdd(&RandNum);
    cout << "We have generated a random number:" << endl; 
    bn_print(&RandNum);
    bool isPrime = rabinmiller(&RandNum, 5);
    cout << "This number ";
    if(isPrime)
        cout << "is ";
    else
        cout << "isn't ";
    cout << "a prime number!" << endl;

    return 0;
}