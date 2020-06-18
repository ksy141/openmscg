#ifndef DELTA_PAIR_BSPLINE_H
#define DELTA_PAIR_BSPLINE_H

#include "delta.h"
#include "bspline.h"

class DeltaPairBSpline : public Delta, public BSpline
{
  public:

    DeltaPairBSpline(class PairList*, int, int, double, double);
    virtual ~DeltaPairBSpline();
    
    virtual void compute();
    
    int type_id;
    class PairList* pair;
};

#endif
