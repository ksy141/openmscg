#ifndef TALBLE_ANGLE_BSPLINE_H
#define TALBLE_ANGLE_BSPLINE_H

#include "table.h"
#include "bspline.h"

class TableAngleBSpline : public Table, public BSpline
{
  public:
    
    TableAngleBSpline(class BondList*, int,
        int, double, double, double);
    virtual ~TableAngleBSpline();

    virtual void compute();
    virtual void dump(double*, double, double, int);
    
    class Topology *top;
    class BondList *lst;
    int type_id;
};

#endif
