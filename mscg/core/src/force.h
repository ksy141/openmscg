#ifndef FORCE_H
#define FORCE_H

#include "pair_list.h"
#include "bond_list.h"

class Force
{
  public:

    static void compute_pair(PairList*, float*, float*);
    static void compute_bond(BondList*, float*, float*);
    static void compute_angle(BondList*, float*, float*);
    static void compute_dihedral(BondList*, float*, float*);
};

#endif