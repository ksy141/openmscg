#ifndef MODEL_ANGLE_BSPLINE_H
#define MODEL_ANGLE_BSPLINE_H

#include "model.h"
#include "bspline.h"

class ModelAngleBSpline : public Model, public BSpline
{
  public:
    
    ModelAngleBSpline(double, double, double, int);
    virtual ~ModelAngleBSpline();

    virtual void compute_rem();
    virtual void compute_fm();
    virtual void get_table(double *params, double* in, double* out, int size);
};

#endif
