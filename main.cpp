/**
#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
using namespace std;


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
#include<iostream>
#include<cstdio>
#include<cstring>
#include<algorithm>
#include<cmath>
using namespace std;
#define ll long long
const ll P=998244353,g=3,G=332748118;
const int N=5e6+5;
ll a[N],b[N];
int order[N];
int p=1,L;
ll qpow(ll a,ll b)
{
    ll ans=1;
    while (b>0)
    {
        if (b&1) ans=ans*a%P;
        b>>=1;
        a=a*a%P;
    }
    return ans;
}
void NTT(ll *a,int inv)
{
	for (int i=0;i<p;i++) 
		if (i<order[i]) swap(a[i],a[order[i]]);
	for (int l=1;l<p;l<<=1)
	{
		ll Wn=qpow(inv==1?g:G,(P-1)/(l<<1));
		for (int R=l<<1,j=0;j<p;j+=R)
		{
			ll w=1;
			for (int k=0;k<l;k++,w=w*Wn%P)
			{
				ll x=a[j+k],y=w*a[j+l+k]%P;
				a[j+k]=(x+y)%P;
				a[j+l+k]=(x-y+P)%P;
			}
		}
	}	
} 
int main()
{
	int n,m;
	scanf("%d%d",&n,&m);
	for (int i=0;i<=n;i++) scanf("%lld",&a[i]),a[i]=(a[i]+P)%P;
	for (int i=0;i<=m;i++) scanf("%lld",&b[i]),b[i]=(b[i]+P)%P;
	while (p<=n+m) p<<=1,L++;
	for (int i=0;i<p;i++)
		order[i]=(order[i>>1]>>1)|((i&1)<<(L-1));
	NTT(a,1),NTT(b,1);
	for (int i=0;i<=p;i++) a[i]=a[i]*b[i]%P;
	NTT(a,-1);
    ll inv=qpow(p,P-2);
	for (int i=0;i<=n+m;i++)
		printf("%d ",a[i]*inv%P);
	return 0;
}
