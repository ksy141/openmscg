#ifndef BSPLINE_H
#define BSPLINE_H

#include "gsl/gsl_bspline.h"
#include "gsl/gsl_blas.h"

class BSpline
{
  public:

    BSpline(int, double, double, double, int scale_flag = 0);
    virtual ~BSpline();

    gsl_bspline_workspace *bw;
    gsl_vector *B;
    
    gsl_vector *B0, *B1, *D0, *D1;
    size_t start0, start1, end0, end1;
    
    int order, ncoeff, scale_flag;
    double xmin, xmax, ddx, resolution;

    double *table_bvalue;
    int *table_istart, *table_nn;

    void setup_cache(double dx_factor = 0.001);
    void eval_coeffs(double, double**, size_t*, int*);
    void eval(double*, double*, double, double, int);
};

#endif
