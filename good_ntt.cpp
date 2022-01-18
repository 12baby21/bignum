#include <iostream>
#include "bn.h"
#include "time.h"
#include "common.h"
#include <limits.h>
#include <cmath>
#include <algorithm>

using namespace std;

#define LL long long

const int mod = 2013265921, g = 31;

struct NTT
{
    void change(int y[], int len) // len必须是2的幂次
    {
        int i, j, k;
        for (i = 1, j = len >> 1; i < len - 1; i++)
        {
            if (i < j)
                swap(y[i], y[j]);
            k = len >> 1;
            while (j >= k)
            {
                j -= k;
                k /= 2;
            }
            if (j < k)
                j += k;
        }       
    }

    void dft(int y[], int len, int on) // len必须是2的幂次
    {
        change(y, len);
        for (int h = 2; h <= len; h <<= 1)
        {
            int wn = qmod(g, (mod - 1) / h, mod);
            if (on == -1)
            {
                wn = qmod(wn, mod - 2, mod);
            }
            for (int j = 0; j < len; j += h)
            {
                int w = 1;
                for (int k = j; k < j + h / 2; k++)
                {
                    LL u = y[k];
                    // cout << "u: "<< u <<endl;
                    LL t = (LL)w * y[k + h / 2] % mod;
                    // cout << "t: "<< t <<endl;
                    y[k] = (u + t) % mod;
                    y[k + h / 2] = (u - t + mod) % mod;
                    w = (LL)w * wn % mod;
                }
            }
        }
        if (on == -1)
        {
            int inv = qmod(len, mod - 2, mod);

            for (int i = 0; i < len; i++)
            {
                y[i] = (LL)y[i] * inv % mod;
            }
        }

    } // on为1时是把系数表示转化为点值表示即DFT，为-1时即逆DFT，还原为系数表示

    inline void mul(int x1[], int n1, int x2[], int n2, int res[], int &n)
    {
        n = 1;
        int _n = n1 + n2 - 1;
        while (n < _n)
            n <<= 1; //把len填充到2的幂次

        dft(x1, n, 1); //转化为点值
        dft(x2, n, 1); //转化为点值
        for (int i = 0; i < n; i++)
        {
            res[i] = (LL)x1[i] * x2[i] % mod;
        }
        dft(res, n, -1); //转化为系数

        for (int i = 0; i < n; i++)
        {
            res[i + 1] += res[i] / 32768;
            res[i] = res[i] % 32768;
        }
        // n=n1+n2-1;
        while (res[n - 1] == 0)
            n--;
    }
} NTT;

const int maxx = (1 << 21) + 5;
int x1[64], x2[64];
int res[maxx];
/**
int main()
{
    int n1 = 32;
    int n2 = 32;
    int n;
    x1[0] = 32767;
    x2[0] = 32767;
    x1[1] = 2;
    x2[1] = 2;


    NTT.mul(x1, n1, x2, n2, res, n);
    for (int i = n - 1; i >= 0; i--)
        printf("%d ", res[i]);
    if (n <= 0)
        cout << 0;
    return 0;
}
**/

