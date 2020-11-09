#include "bond_list.h"

#include <cmath>
#include <cstdio>
#include <cassert>

BondList::BondList(
    int nbonds, int *bond_types, vec2i *bonds, 
    int nangles, int *angle_types, vec3i *angles,
    int ndihedrals, int *dihedral_types, vec4i *dihedrals)
{
    this->nbonds = nbonds;
    this->bond_types = bond_types;
    this->bond_atoms = bonds;
    
    this->nangles = nangles;
    this->angle_types = angle_types;
    this->angle_atoms = angles;
    
    this->ndihedrals = ndihedrals;
    this->dihedral_types = dihedral_types;
    this->dihedral_atoms = dihedrals;
    
    #define _ptr(p) if(nbonds>0) p##_bond = new float[nbonds]; else p##_bond = 0
    _ptr(dx); _ptr(dy); _ptr(dz); _ptr(dr);
    #undef _ptr
    
    #define _ptr(p) if(nangles>0) p##_angle = new float[nangles]; else p##_angle = 0
    _ptr(dx1); _ptr(dy1); _ptr(dz1); 
    _ptr(dx2); _ptr(dy2); _ptr(dz2);
    _ptr(a11); _ptr(a12); _ptr(a22);
    _ptr(theta);
    #undef _ptr
    
    #define _ptr(p) if(ndihedrals>0) p##_dihedral = new float[ndihedrals]; else p##_dihedral = 0
    _ptr(dpd1x); _ptr(dpd1y); _ptr(dpd1z);
    _ptr(dpd2x); _ptr(dpd2y); _ptr(dpd2z);
    _ptr(dpd3x); _ptr(dpd3y); _ptr(dpd3z);
    _ptr(dpd4x); _ptr(dpd4y); _ptr(dpd4z);
    _ptr(phi);
    #undef _ptr
}

BondList::~BondList()
{
    #define _ptr(p) if(p##_bond) delete [] p##_bond
    _ptr(dx); _ptr(dy); _ptr(dz); _ptr(dr);
    #undef _ptr
    
    #define _ptr(p) if(p##_angle) delete [] p##_angle
    _ptr(dx1); _ptr(dy1); _ptr(dz1); 
    _ptr(dx2); _ptr(dy2); _ptr(dz2);
    _ptr(a11); _ptr(a12); _ptr(a22);
    _ptr(theta);
    #undef _ptr
    
    #define _ptr(p) if(p##_dihedral) delete [] p##_dihedral
    _ptr(dpd1x); _ptr(dpd1y); _ptr(dpd1z);
    _ptr(dpd2x); _ptr(dpd2y); _ptr(dpd2z);
    _ptr(dpd3x); _ptr(dpd3y); _ptr(dpd3z);
    _ptr(dpd4x); _ptr(dpd4y); _ptr(dpd4z);
    _ptr(phi);
    #undef _ptr
}

void BondList::build(vec3f box, vec3f *x)
{
    if(bond_atoms) build_bonds(box, x);
    if(angle_atoms) build_angles(box, x); 
    if(dihedral_atoms) build_dihedrals(box, x);
}

void BondList::build_bonds(vec3f box, vec3f *x)
{    
    for(int i=0; i<nbonds; i++)
    {
        int ia = bond_atoms[i][0], ib = bond_atoms[i][1];
        float r[3];
        vector_sub(r, x[ib], x[ia]);
        
        float r2 = vector_dot(r, r);
        float dr = sqrt(r2);
        
        dx_bond[i] = r[0];
        dy_bond[i] = r[1];
        dz_bond[i] = r[2];
        dr_bond[i] = dr;
    }
}

void BondList::build_angles(vec3f box, vec3f *x)
{
    for(int i=0; i<nangles; i++)
    {
        int ia = angle_atoms[i][0], ib = angle_atoms[i][1], ic = angle_atoms[i][2];
        
        float d1[3], d2[3];
        vector_sub(d1, x[ia], x[ib]);
        vector_sub(d2, x[ic], x[ib]);        
        
        float rsq1 = vector_dot(d1, d1);
        float dr1 = sqrt(rsq1);
        
        float rsq2 = vector_dot(d2, d2);
        float dr2 = sqrt(rsq2);
        
        float c = vector_dot(d1, d2);
        c /= dr1*dr2;
        
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        
        float s = sqrt(1.0 - c*c);
        if (s < 0.0001) s = 0.0001;
        s = 1.0/s;
        
        theta_angle[i] = acos(c);
        dx1_angle[i] = d1[0]; dy1_angle[i] = d1[1]; dz1_angle[i] = d1[2];
        dx2_angle[i] = d2[0]; dy2_angle[i] = d2[1]; dz2_angle[i] = d2[2];
        a11_angle[i] = s * c / rsq1;
        a12_angle[i] = -s / (dr1 * dr2);
        a22_angle[i] = s * c / rsq2;
    }
}

void BondList::build_dihedrals(vec3f box, vec3f *x)
{
    for(int i=0; i<ndihedrals; i++)
    {
        int i1 = dihedral_atoms[i][0];
        int i2 = dihedral_atoms[i][1];
        int i3 = dihedral_atoms[i][2];
        int i4 = dihedral_atoms[i][3];
        
        float b1[3], b2[3], b3[3];
        vector_sub(b1, x[i2], x[i1]);
        vector_sub(b2, x[i3], x[i2]);
        vector_sub(b3, x[i4], x[i3]);
        
        float n1[3], n2[3];
        vector_prod(n1, b2, b1);
        vector_prod(n2, b2, b3);
        vector_norm(n1);
        vector_norm(n2);
                
        float c = -1.0 * vector_dot(n1, n2);
        if(c > 1.0) c = 1.0; else if(c < -1.0) c = -1.0;
        
        float phi = acos(c);
        if(vector_dot(n1, b3) > 0.0) phi = M_PI * 2.0 - phi;
        
        phi_dihedral[i] = phi;
        
        float p12 = vector_dot(b1, b2);
        float p23 = vector_dot(b2, b3);
        
        float L2sq = vector_dot(b2, b2);
        float L2 = sqrt(L2sq);
        float L2sqinv = L2sq>0.0?1.0/L2sq:0.0;
        float L2inv = sqrt(L2sqinv);
                
        float p1to2[3], p3to2[3], p1to2r[3], p3to2r[3];
        vector_scale(p1to2, b2, p12 * L2inv);
        vector_scale(p3to2, b2, p23 * L2inv);
        vector_sub(p1to2r, b1, p1to2);
        vector_sub(p3to2r, b3, p3to2);
        
        float L12 = sqrt(vector_dot(p1to2r, p1to2r));
        float L12inv = L12>0.0?1.0/L12:0.0;
        
        float L34 = sqrt(vector_dot(p3to2r, p3to2r));
        float L34inv = L34>0.0?1.0/L34:0.0;
        
        dpd1x_dihedral[i] = n1[0] * L12inv;
        dpd1y_dihedral[i] = n1[1] * L12inv;
        dpd1z_dihedral[i] = n1[2] * L12inv;
        
        dpd4x_dihedral[i] = n2[0] * L34inv;
        dpd4y_dihedral[i] = n2[1] * L34inv;
        dpd4z_dihedral[i] = n2[2] * L34inv;
        
        float dp123dx2 = -L2inv * (L2 + p12 * L2inv);
        float dp234dx2 = L2inv * p23 * L2inv;
        
        float dp123dx3 = L2inv * p12 * L2inv;
        float dp234dx3 = -L2inv * (L2 + p23 * L2inv);
        
        dpd2x_dihedral[i] = dp123dx2 * dpd1x_dihedral[i] + dp234dx2 * dpd4x_dihedral[i];
        dpd2y_dihedral[i] = dp123dx2 * dpd1y_dihedral[i] + dp234dx2 * dpd4y_dihedral[i];
        dpd2z_dihedral[i] = dp123dx2 * dpd1z_dihedral[i] + dp234dx2 * dpd4z_dihedral[i];
        
        dpd3x_dihedral[i] = dp123dx3 * dpd1x_dihedral[i] + dp234dx3 * dpd4x_dihedral[i];
        dpd3y_dihedral[i] = dp123dx3 * dpd1y_dihedral[i] + dp234dx3 * dpd4y_dihedral[i];
        dpd3z_dihedral[i] = dp123dx3 * dpd1z_dihedral[i] + dp234dx3 * dpd4z_dihedral[i];
    }
}
