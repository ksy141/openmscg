#ifndef TALBLE_BOND_BSPLINE_H
#define TALBLE_BOND_BSPLINE_H

#include "table.h"
#include "bspline.h"

class TableBondBSpline : public Table, public BSpline
{
  public:
    
    TableBondBSpline(class BondList*, int,
        int, double, double, double);
    virtual ~TableBondBSpline();

    virtual void compute();
    virtual void dump(double*, double, double, int);
    
    class Topology *top;
    class BondList *lst;
    int type_id;
};

#endif
