#ifndef MODEL_PAIR_BSPLINE_H
#define MODEL_PAIR_BSPLINE_H

#include "model.h"
#include "bspline.h"

class ModelPairBSpline : public Model, public BSpline
{
  public:

    ModelPairBSpline(class PairList*, int, int, double, double);
    virtual ~ModelPairBSpline();
    
    virtual void compute_dudl();
    virtual void compute_dfdl();
    virtual void compute_etable(double *params, double* in, double* out, int size);
    virtual void compute_ftable(double *params, double* in, double* out, int size);
};

#endif
