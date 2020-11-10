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
        if(vector_dot(n1, b3) > 0.0) phi = - phi;
        
        phi_dihedral[i] = phi;
        
        // derivatives
        // ref1: https://salilab.org/modeller/9v6/manual/node436.html
        // ref2: https://grigoryanlab.org/docs/dynamics_derivatives.pdf
        
        float r_ij[3], r_kj[3], r_kl[3];
        float r_mj[3], r_nk[3];
        float f1[3], f2[3], f3[3], f4[3];
        float f2a[3], f2b[3], f3a[3], f3b[3];
        
        vector_sub(r_ij, x[i2], x[i1]);
        vector_sub(r_kj, x[i2], x[i3]);
        vector_sub(r_kl, x[i4], x[i3]);
        
        vector_prod(r_mj, r_ij, r_kj);
        vector_prod(r_nk, r_kj, r_kl);
        
        float l2_kj = vector_dot(r_kj, r_kj);
        float l_kj = sqrt(l2_kj);
        
        float s1 = l_kj / vector_dot(r_mj, r_mj);
        vector_scale(f1, r_mj, s1);
        
        float s4 = - l_kj / vector_dot(r_nk, r_nk);
        vector_scale(f4, r_nk, s4);
        
        float s2a = vector_dot(r_ij, r_kj) / l2_kj - 1.0;
        vector_scale(f2a, f1, s2a);
        float s2b = vector_dot(r_kl, r_kj) / l2_kj * (-1.0);
        vector_scale(f2b, f4, s2b);
        vector_add(f2, f2a, f2b);
        
        float s3a = vector_dot(r_kl, r_kj) / l2_kj - 1.0;
        vector_scale(f3a, f4, s3a);
        float s3b = vector_dot(r_ij, r_kj) / l2_kj * (-1.0);
        vector_scale(f3b, f1, s3b);
        vector_add(f3, f3a, f3b);
        
        dpd1x_dihedral[i] = f1[0];
        dpd1y_dihedral[i] = f1[1];
        dpd1z_dihedral[i] = f1[2];
        
        dpd2x_dihedral[i] = f2[0];
        dpd2y_dihedral[i] = f2[1];
        dpd2z_dihedral[i] = f2[2];
        
        dpd3x_dihedral[i] = f3[0];
        dpd3y_dihedral[i] = f3[1];
        dpd3z_dihedral[i] = f3[2];
        
        dpd4x_dihedral[i] = f4[0];
        dpd4y_dihedral[i] = f4[1];
        dpd4z_dihedral[i] = f4[2];
    }
}
