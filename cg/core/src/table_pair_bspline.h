#ifndef TALBLE_PAIR_BSPLINE_H
#define TALBLE_PAIR_BSPLINE_H

#include "table.h"
#include "bspline.h"

class TablePairBSpline : public Table, public BSpline
{
  public:

    TablePairBSpline(class PairList*, int, int, double, double);
    virtual ~TablePairBSpline();
    
    virtual void compute();
    virtual void dump(double*, double, double, int);
    
    int type_id;
    class PairList* pair;
};

#endif
