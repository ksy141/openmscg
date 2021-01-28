#include "tables.h"
#include <cmath>
#include <cstdio>

Tables::Tables(int maxtable)
{
    tbls = new Table[maxtable];
    for(int i=0; i<maxtable; i++) tbls[i].n = 0;
}

Tables::~Tables()
{
    delete [] tbls;
}

void Tables::compute(long n, int *tids, float *x, float *e, float *f)
{
    for(long i=0; i<n; i++) 
    {
        int tid = tids[i];
        if(tbls[tid].n == 0) continue;
        
        int it = round((x[i] - tbls[tid].min) / tbls[tid].inc);
        if(it<0) it = 0; else if(it>=tbls[tid].n) it = tbls[tid].n;
        
        e[i] += tbls[tid].efac[it];
        f[i] += tbls[tid].ffac[it];
    }
}