#ifndef MODEL_NB3B_SW_H
#define MODEL_NB3B_SW_H

#include "model.h"

class ModelNB3BSW : public Model
{
  public:

    ModelNB3BSW(double, double, double, double, double);
    virtual ~ModelNB3BSW();

    virtual void compute_rem();
    virtual void compute_fm();
    virtual void get_table(double *params, double* in, double* out, int size);
    void setup_ex(int, int);
    void compute_fm_one(int, int, int, float*, float*, float, float);

    int tid_i, tid_ij, tid_ik;
    float gamma_ij, gamma_ik, a_ij, a_ik, cos0;
};

#endif
