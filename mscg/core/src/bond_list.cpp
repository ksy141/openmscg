#include "bond_list.h"

#include <cmath>
#include <cstdio>
#include <cassert>

BondList::BondList(int nbonds, vec2i *bonds, int nangles, vec3i *angles, int ndihedrals, vec4i *dihedrals)
{
    this->nbonds = nbonds;
    this->bond_atoms = bonds;
    
    this->nangles = nangles;
    this->angle_atoms = angles;
    
    this->ndihedrals = ndihedrals;
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
        vector_sub(d, x[ib], x[ia]);
        
        float r2 = dx*dx + dy*dy + dz*dz;
        float dr = sqrt(r2);
        
        dx_bond[i] = dx;
        dy_bond[i] = dy;
        dz_bond[i] = dz;
        dr_bond[i] = dr;
    }
}

void BondList::build_angles(vec3f box, vec3f *x)
{
    for(int i=0; i<nangles; i++)
    {
        int ia = angle_atoms[i][0], ib = angle_atoms[i][1], ic = angle_atoms[i][2];
        
        vector_sub(d1, x[ia], x[ib]);
        vector_sub(d2, x[ic], x[ib]);        
        
        float rsq1 = d1x*d1x + d1y*d1y + d1z*d1z;
        float dr1 = sqrt(rsq1);
        
        float rsq2 = d2x*d2x + d2y*d2y + d2z*d2z;
        float dr2 = sqrt(rsq2);
        
        float c = d1x*d2x + d1y*d2y + d1z*d2z;
        c /= dr1*dr2;
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        
        float s = sqrt(1.0 - c*c);
        if (s < 0.0001) s = 0.0001;
        s = 1.0/s;
        
        theta_angle[i] = acos(c);
        dx1_angle[i] = d1x; dy1_angle[i] = d1y; dz1_angle[i] = d1z;
        dx2_angle[i] = d2x; dy2_angle[i] = d2y; dz2_angle[i] = d2z;
        a11_angle[i] = s * c / rsq1;
        a12_angle[i] = -s / (dr1 * dr2);
        a22_angle[i] = s * c / rsq2;
    }
}

void BondList::build_dihedrals(vec3f box, vec3f *x)
{
    
}