#ifndef NB3B_H
#define NB3B_H

#include "pair_list.h"

class NB3B {};

class NB3B_SW : public NB3B
{
public:
    NB3B_SW(int, int, int, float, float, float, float, float, float);
    virtual ~NB3B_SW();

    void compute(PairList*, float*, float*, float*);
    float compute_one(float, float, float*, float*, float*, float*, float*);
    
    int type_i, type_j, type_k;
    float lambda, cos0;
    float gamma_ij, gamma_ik, a_ij, a_ik;
};

#endif
