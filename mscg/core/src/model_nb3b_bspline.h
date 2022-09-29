#ifndef MODEL_NB3B_BSPLINE_H
#define MODEL_NB3B_BSPLINE_H

#include "model.h"
#include "bspline.h"

class ModelNB3BBSpline : public Model, public BSpline
{
  public:

    ModelNB3BBSpline(double, double, double, int);
    virtual ~ModelNB3BBSpline();

    virtual void compute_rem();
    virtual void compute_fm();
    virtual void get_table(double *params, double* in, double* out, int size);

    void setup_ex(int, double, double, int, double, double);
    void compute_fm_one(int, int, int, float*, float*, float, float);

    int tid_i, tid_ij, tid_ik;
    float gamma_ij, gamma_ik, a_ij, a_ik;
    bool using_theta;
};

#endif
