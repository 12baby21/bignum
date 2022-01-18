#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
#include <limits.h>
#include<cmath>
#include<algorithm>

using namespace std;


#define ll long long



const int mod = 2013265921, g = 31;

int P(int a, int b)
{
    int ans = 1;
    while(b)
    {
        if(b & 1)
            ans = (ll)ans * a % mod;
        a = (ll)a * a % mod;
        b >>= 1;
    }
    return ans;
}

struct NTT
{
    void change(int y[], int len)//len必须是2的幂次
    {
        int i, j, k;
        for(i = 1, j = len>>1; i < len-1; i++)
        {
            if(i < j)
                swap(y[i], y[j]);
            k = len>>1;
            while(j >= k)
            {
                j -= k;
                k /= 2;
            }
            if(j < k)
                j += k;
        }
    }

    void dft(int y[], int len, int on)//len必须是2的幂次
    {
        change(y, len);
        for(int h=2; h <= len; h<<=1)
        {
            int wn = P(g, on==1? (mod-1)/h: mod-1-(mod-1)/h);
            for(int j = 0; j < len; j += h)
            {
                int w = 1;
                for(int k = j; k < j+h/2; k++)
                {
                    int u = y[k];
                    int t = (ll)w * y[k+h/2] % mod;
                    y[k] = (u+t)%mod;
                    y[k+h/2] = (u-t+mod) % mod;
                    w = (ll)w*wn%mod;
                }
            }
        }
        if(on==-1)
        {
            int inv = P(len, mod-2);
            for(int i = 0 ;i < len; i++)
                y[i] = (ll)y[i] * inv % mod;
        }

    }//on为1时是把系数表示转化为点值表示即DFT，为-1时即逆DFT，还原为系数表示

    inline void mul(int x1[], int n1, int x2[], int n2, int res[], int& n)
    {
        n=1;
        int _n=n1+n2-1;
        while(n < _n)
            n <<= 1;//把len填充到2的幂次
        cout<<n<<endl;
        dft(x1,n,1);//转化为点值
        dft(x2,n,1);//转化为点值
        for(int i = 0; i < n; i++)
            res[i]=(ll)x1[i] * x2[i] % mod;
        dft(res,n,-1);//转化为系数
        for(int i = 0; i < n; i++)
        {
            res[i+1] += res[i] / 32768;
            res[i] = res[i] % 32768;
        }
        //n=n1+n2-1;
        while(res[n-1]==0)
            n--;
    }
}NTT;

const int maxx=(1<<21)+5;
int x1[32], x2[32];
int res[maxx];

int main()
{
    int n1 = 16;
    int n2 = 16;
    int n;
    x1[0] = 32760;
    x2[0] = 32760;

    NTT.mul(x1,n1,x2,n2,res,n);
    for(int i = n-1; i >= 0; i--)
        printf("%d ", res[i]);
    if(n <= 0 )
        cout<< 0;
    return 0;
}



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
