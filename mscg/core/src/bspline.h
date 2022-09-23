#ifndef BSPLINE_H
#define BSPLINE_H

#include "gsl/gsl_bspline.h"
#include "gsl/gsl_blas.h"

#define SMALL 0.000001

class BSpline
{
  public:

    BSpline(int, double, double, double);
    virtual ~BSpline();

    gsl_bspline_workspace *bw;
    gsl_vector *B, *D;
    gsl_matrix *m, *mt;

    gsl_vector *B0, *B1, *D0, *D1;
    size_t start0, start1, end0, end1;

    int order, nbreak, ncoeff;
    double xmin, xmax, ddx, resolution;

    double *table_bvalue;
    int *table_istart, *table_nn;

    static inline int get_nbreak(double xmin, double xmax, double resolution)
    {
        return static_cast<int>(ceil((xmax - xmin)/resolution)) + 1;
    }

    void setup_cache(double dx_factor = 0.001);
    void eval_coeffs(double, double**, size_t*, int*);
    void eval_derivs(double, double**, size_t*, int*);
    void eval(double*, double*, double, double, int);
    void eval(double*, double*, double*, int);

    enum ExtrapolationMethod {LINEAR, TRUNC, CAP};
    ExtrapolationMethod extrapolation;
};

#endif
