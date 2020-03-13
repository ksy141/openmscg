#include "bond_list.h"

#include <cmath>
#include <cstdio>

BondList::BondList(Topology *top, Traj *trj)
{
    this->top = top;
    this->trj = trj;
    
    dx_bond = dy_bond = dz_bond = dr_bond = 0;
    
    #define _ptr(p) if(top->nbonds>0) p##_bond = new float[top->nbonds]; else p##_bond = 0
    _ptr(dx); _ptr(dy); _ptr(dz); _ptr(dr);
    #undef _ptr
    
    #define _ptr(p) if(top->nangls>0) p##_angl = new float[top->nangls]; else p##_angl = 0
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
    
    #define _ptr(p) if(p##_angl) delete [] p##_angl
    _ptr(dx1); _ptr(dy1); _ptr(dz1); 
    _ptr(dx2); _ptr(dy2); _ptr(dz2);
    _ptr(a11); _ptr(a12); _ptr(a22);
    _ptr(theta);
    #undef _ptr
}

void BondList::build()
{
    for(int d=0; d<3; d++) 
    {
        box[d] = trj->box[d];
        hbox[d] = 0.5 * box[d];
    }
    
    build_bonds();
    build_angls(); 
    build_dihes();
}

void BondList::build_bonds()
{    
    int *atom1 = top->bond_atom1;
    int *atom2 = top->bond_atom2;
    Vec *x = trj->x;

    for(int i=0; i<top->nbonds; i++)
    {
        int ia = atom1[i], ib = atom2[i];
        float dx = x[ib][0] - x[ia][0];
        float dy = x[ib][1] - x[ia][1];
        float dz = x[ib][2] - x[ia][2];
        
        if(dx>hbox[0]) dx-=box[0]; else if(dx<-hbox[0]) dx+=box[0];
        if(dy>hbox[1]) dy-=box[1]; else if(dy<-hbox[1]) dy+=box[1];
        if(dz>hbox[2]) dz-=box[2]; else if(dz<-hbox[2]) dz+=box[2];
        
        float r2 = dx*dx + dy*dy + dz*dz;
        float dr = sqrt(r2);
        
        dx_bond[i] = dx;
        dy_bond[i] = dy;
        dz_bond[i] = dz;
        dr_bond[i] = dr;
    }
}

void BondList::build_angls()
{
    int *atom1 = top->angl_atom1;
    int *atom2 = top->angl_atom2;
    int *atom3 = top->angl_atom3;
    Vec *x = trj->x;

    for(int i=0; i<top->nangls; i++)
    {
        int ia = atom1[i], ib = atom2[i], ic = atom3[i];
        
        float dx1 = x[ia][0] - x[ib][0];
        float dy1 = x[ia][1] - x[ib][1];
        float dz1 = x[ia][2] - x[ib][2];
        
        if(dx1>hbox[0]) dx1-=box[0]; else if(dx1<-hbox[0]) dx1+=box[0];
        if(dy1>hbox[1]) dy1-=box[1]; else if(dy1<-hbox[1]) dy1+=box[1];
        if(dz1>hbox[2]) dz1-=box[2]; else if(dz1<-hbox[2]) dz1+=box[2];
        
        float dx2 = x[ic][0] - x[ib][0];
        float dy2 = x[ic][1] - x[ib][1];
        float dz2 = x[ic][2] - x[ib][2];
        
        if(dx2>hbox[0]) dx2-=box[0]; else if(dx2<-hbox[0]) dx2+=box[0];
        if(dy2>hbox[1]) dy2-=box[1]; else if(dy2<-hbox[1]) dy2+=box[1];
        if(dz2>hbox[2]) dz2-=box[2]; else if(dz2<-hbox[2]) dz2+=box[2];
        
        float rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
        float dr1 = sqrt(rsq1);
        
        float rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
        float dr2 = sqrt(rsq2);
        
        float c = dx1*dx2 + dy1*dy2 + dz1*dz2;
        c /= dr1*dr2;
        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;
        
        float s = sqrt(1.0 - c*c);
        if (s < 0.0001) s = 0.0001;
        s = 1.0/s;
        
        theta_angl[i] = acos(c);
        dx1_angl[i] = dx1; dy1_angl[i] = dy1; dz1_angl[i] = dz1;
        dx2_angl[i] = dx2; dy2_angl[i] = dy2; dz2_angl[i] = dz2;
        a11_angl[i] = s * c / rsq1;
        a12_angl[i] = -s / (dr1 * dr2);
        a22_angl[i] = s * c / rsq2;
    }
}

void BondList::build_dihes()
{

}