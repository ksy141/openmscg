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
    
    float *dx_bond, *dy_bond, *dz_bond, *dr_bond;
    
    float *theta_angl;
    float *dx1_angl, *dy1_angl, *dz1_angl;
    float *dx2_angl, *dy2_angl, *dz2_angl;
    float *a11_angl, *a12_angl, *a22_angl;
    
    // functions
    
    BondList(Topology*);
    virtual ~BondList();
    
    void build(Traj*);
    void build_bonds(Traj*);
    void build_angls(Traj*);
    void build_dihes(Traj*);
};

#endif
