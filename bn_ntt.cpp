#include <bits/stdc++.h>
#include <bn.h>

using namespace std;
typedef long long LL;
const int N = 1 << 18;
const int G = 3, P = (479 << 21) + 1; // G为原根，P为大素数,可以处理2^21次范围

LL quick(LL x, LL n)
{
    LL ret = 1;
    for (; n; n >>= 1)
    {
        if (n & 1)
            ret = ret * x % P;
        x = x * x % P;
    }
    return ret;
}

LL A[N], B[N];
void rader(LL *y, int len)
{
    for (int i = 1, j = len / 2; i < len - 1; i++)
    {
        if (i < j)
            swap(y[i], y[j]);
        int k = len / 2;
        while (j >= k)
        {
            j -= k;
            k /= 2;
        }
        if (j < k)
            j += k;
    }
}
void ntt(bn_ptr y, int len = 32, int op = 1)
{
    int l = 0;
    while ((l << 1) < len)
        ++l;

    for (int r = 1; r <= l; ++r)
    {
        int mid = 1 << (l - r);
        bn t;
        bn g(G);
        bn M(P);
        bn_sub(&M, &M, &one);
        bn_qpow(&t, &g, )
        bn wn;
        
        bn
            bn_qpow(&wn, G, P - 1 / h);
        quick(G, (P - 1) / h);
        if (op == -1)
            wn = quick(wn, P - 2);
        for (int j = 0; j < len; j += h)
        {
            LL w = 1;
            for (int k = j; k < j + h / 2; k++)
            {
                LL u = y[k];
                LL t = w * y[k + h / 2] % P;
                y[k] = (u + t) % P;
                y[k + h / 2] = (u - t + P) % P;
                w = w * wn % P;
            }
        }
    }
    if (op == -1)
    {
        LL inv = quick(len, P - 2);
        for (int i = 0; i < len; i++)
            y[i] = y[i] * inv % P;
    }
}
char s1[100010], s2[100010];
LL ans[N];

int main()
{
    while (scanf("%s%s", s1, s2) != EOF)
    {
        int len1 = strlen(s1);
        int len2 = strlen(s2);
        memset(A, 0, sizeof(A));
        memset(B, 0, sizeof(B));
        for (int i = len1 - 1; i >= 0; i--)
            A[i] = s1[len1 - i - 1] - '0';
        for (int i = len2 - 1; i >= 0; i--)
            B[i] = s2[len2 - i - 1] - '0';
        int len = 1;
        while (len < len1 * 2 || len < len2 * 2)
            len <<= 1;
        ntt(A, len, 1);
        ntt(B, len, 1);
        for (int i = 0; i < len; i++)
            A[i] = A[i] * B[i] % P;
        ntt(A, len, -1);
        memset(ans, 0, sizeof(ans));
        for (int i = 0; i < len; i++)
        {
            ans[i] += A[i];
            if (ans[i] >= 10)
            {
                ans[i + 1] += ans[i] / 10;
                ans[i] %= 10;
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
        cout << s1 << " * " << s2 << " = ";
        for (int i = pos; i >= 0; i--)
        {
            cout << ans[i];
        }
        cout << endl;
    }
    return 0;
}
