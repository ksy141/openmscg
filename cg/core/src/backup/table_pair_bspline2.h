#ifndef TALBLE_PAIR_BSPLINE_H
#define TALBLE_PAIR_BSPLINE_H

#include "table.h"
#include "gsl/gsl_bspline.h"
#include "gsl/gsl_blas.h"

class TablePairBSpline : public Table
{
  public:
    
    gsl_bspline_workspace *bw;
    gsl_vector *B;
    
    TablePairBSpline(class PairList*, int, int, double, double);
    virtual ~TablePairBSpline();
    
    virtual void compute();
    virtual void dump(double*, double, double, int);
    
    int type_id;
    double xmin, xmax, ddx, resolution;
    class PairList* pair;
    int order;
    
    double *table_bvalue;
    int *table_istart, *table_nn;
    
    void setup_cache(double dx_factor = 0.001);
};

#endif
