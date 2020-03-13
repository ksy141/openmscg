#ifndef TALBLE_BOND_BSPLINE_H
#define TALBLE_BOND_BSPLINE_H

#include "table.h"
#include "gsl/gsl_bspline.h"
#include "gsl/gsl_blas.h"

class TableBondBSpline : public Table
{
  public:

    gsl_bspline_workspace *bw;
    gsl_vector *B;

    TableBondBSpline(class BondList*, int,
        int, double, double, double);
    virtual ~TableBondBSpline();

    virtual void compute();
    virtual void dump(double*, double, double, int);

    void setup_cache(double dx_factor=0.001);

    class Topology *top;
    class BondList *lst;
    int type_id, order;
    double resolution, xmin, xmax;

    double *table_bvalue;
    int *table_istart, *table_nn;
    double ddx;
};

#endif
