#ifndef BOND_LIST_H
#define BOND_LIST_H

#include "defs.h"

class BondList
{
  public:
    
    // bonds
    int nbonds;
    vec2i *bond_atoms;
    int *bond_types;
    float *dr_bond;
    float *dx_bond, *dy_bond, *dz_bond;
    
    // angles
    int nangles;
    vec3i *angle_atoms;
    int *angle_types;
    float *theta_angle;
    float *dx1_angle, *dy1_angle, *dz1_angle;
    float *dx2_angle, *dy2_angle, *dz2_angle;
    float *a11_angle, *a12_angle, *a22_angle;
    
    // dihedrals
    int ndihedrals;
    vec4i *dihedral_atoms;
    int *dihedral_types;
    float *phi_dihedral;
    
    float *dpd1x_dihedral, *dpd1y_dihedral, *dpd1z_dihedral;
    float *dpd2x_dihedral, *dpd2y_dihedral, *dpd2z_dihedral;
    float *dpd3x_dihedral, *dpd3y_dihedral, *dpd3z_dihedral;
    float *dpd4x_dihedral, *dpd4y_dihedral, *dpd4z_dihedral;
    
    // functions
    
    BondList(int nbonds, int*, vec2i *bonds, int nangles, int*, vec3i *angles, int ndihedrals, int*, vec4i *dihedrals);
    virtual ~BondList();
    
    void build(vec3f box, vec3f *x);
    void build_bonds(vec3f box, vec3f *x);
    void build_angles(vec3f box, vec3f *x);
    void build_dihedrals(vec3f box, vec3f *x);
};

#endif
