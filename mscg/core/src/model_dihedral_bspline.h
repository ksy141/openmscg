#ifndef MODEL_DIHEDRAL_BSPLINE_H
#define MODEL_DIHEDRAL_BSPLINE_H

#include "model.h"
#include "bspline.h"

class ModelDihedralBSpline : public Model, public BSpline
{
  public:
    
    ModelDihedralBSpline(double, double, double, int);
    virtual ~ModelDihedralBSpline();

    virtual void compute_rem();
    virtual void compute_fm();
    virtual void get_table(double *params, double* in, double* out, int size);
};

#endif
