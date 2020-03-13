#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "traj.h"
#include "topology.h"

class BondList
{
  public:
    
    // settings
    
    Vec box, hbox;
    
    Topology* top;
    Traj* trj;
    
    float *dx_bond, *dy_bond, *dz_bond, *dr_bond;
    
    float *theta_angl;
    float *dx1_angl, *dy1_angl, *dz1_angl;
    float *dx2_angl, *dy2_angl, *dz2_angl;
    float *a11_angl, *a12_angl, *a22_angl;
    
    // functions
    
    BondList(Topology*, Traj*);
    virtual ~BondList();
    
    void build();
    void build_bonds();
    void build_angls();
    void build_dihes();
};

#endif
