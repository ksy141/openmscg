#ifndef MODEL_PAIR_BSPLINE_H
#define MODEL_PAIR_BSPLINE_H

#include "model.h"
#include "bspline.h"

class ModelBondBSpline : public Model, public BSpline
{
  public:

    ModelBondBSpline(double, double, double, int);
    virtual ~ModelBondBSpline();
    
    virtual void compute_rem();
    virtual void compute_fm();
    virtual void get_table(double *params, double* in, double* out, int size);
};

#endif
