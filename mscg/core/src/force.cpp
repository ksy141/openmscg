#include "force.h"

void Force::compute_pair(PairList *pair, float *dU, float *f)
{
    int* ilist = pair->ilist;
    int* jlist = pair->jlist;
    int* tlist = pair->tlist;
    
    float* dx = pair->dxlist;
    float* dy = pair->dylist;
    float* dz = pair->dzlist;
    float* dr = pair->drlist;

    for(int p=0; p<pair->npairs; p++)
    {
        float ff = - dU[p] / dr[p];
        float fx = ff * dx[p];
        float fy = ff * dy[p];
        float fz = ff * dz[p];
        
        float *fi = f + ilist[p] * 3;
        fi[0] += fx;
        fi[1] += fy;
        fi[2] += fz;
        
        float *fj = f + jlist[p] * 3;
        fj[0] -= fx;
        fj[1] -= fy;
        fj[2] -= fz;
    }
}

void Force::compute_bond(BondList *blist, float *dU, float *f)
{
    vec2i *atoms = blist->bond_atoms;
    float *dx = blist->dx_bond;
    float *dy = blist->dy_bond;
    float *dz = blist->dz_bond;
    float *dr = blist->dr_bond;

    for(int i=0; i<blist->nbonds; i++)
    {
        float ff = - dU[i] / dr[i];
        float fx = ff * dx[i];
        float fy = ff * dy[i];
        float fz = ff * dz[i];
                
        float *fi = f + atoms[i][0] * 3;
        fi[0] += fx;
        fi[1] += fy;
        fi[2] += fz;
        
        float *fj = f + atoms[i][1] * 3;
        fj[0] -= fx;
        fj[1] -= fy;
        fj[2] -= fz;
    }
}

void Force::compute_angle(BondList *lst, float *dU, float *f)
{
    vec3i *atoms = lst->angle_atoms;
    float *dx1 = lst->dx1_angle;
    float *dy1 = lst->dy1_angle;
    float *dz1 = lst->dz1_angle;
    float *dx2 = lst->dx2_angle;
    float *dy2 = lst->dy2_angle;
    float *dz2 = lst->dz2_angle;
    float *a11 = lst->a11_angle;
    float *a12 = lst->a12_angle;
    float *a22 = lst->a22_angle;

    for(int i=0; i<lst->nangles; i++)
    {
        float ff = dU[i];
        float f1x = a11[i]*dx1[i] + a12[i]*dx2[i];
        float f1y = a11[i]*dy1[i] + a12[i]*dy2[i];
        float f1z = a11[i]*dz1[i] + a12[i]*dz2[i];
        float f3x = a22[i]*dx2[i] + a12[i]*dx1[i];
        float f3y = a22[i]*dy2[i] + a12[i]*dy1[i];
        float f3z = a22[i]*dz2[i] + a12[i]*dz1[i];
        
        float *f1 = f + atoms[i][0] * 3;
        f1[0] += f1x * ff;
        f1[1] += f1y * ff;
        f1[2] += f1z * ff;
        
        float *f3 = f + atoms[i][2] * 3;
        f3[0] += f3x * ff;
        f3[1] += f3y * ff;
        f3[2] += f3z * ff;
        
        float *f2 = f + atoms[i][1] * 3;
        f2[0] -= (f1x + f3x) * ff;
        f2[1] -= (f1y + f3y) * ff;
        f2[2] -= (f1z + f3z) * ff;
    }
}

void Force::compute_dihedral(BondList *lst, float *dU, float *f)
{
    vec4i *atoms = lst->dihedral_atoms;
    float *dpd1x = lst->dpd1x_dihedral;
    float *dpd1y = lst->dpd1y_dihedral;
    float *dpd1z = lst->dpd1z_dihedral;
    float *dpd2x = lst->dpd2x_dihedral;
    float *dpd2y = lst->dpd2y_dihedral;
    float *dpd2z = lst->dpd2z_dihedral;
    float *dpd3x = lst->dpd3x_dihedral;
    float *dpd3y = lst->dpd3y_dihedral;
    float *dpd3z = lst->dpd3z_dihedral;
    float *dpd4x = lst->dpd4x_dihedral;
    float *dpd4y = lst->dpd4y_dihedral;
    float *dpd4z = lst->dpd4z_dihedral;

    for(int i=0; i<lst->ndihedrals; i++)
    {
        float ff = dU[i];
        
        float* f1 = f + atoms[i][0] * 3;        
        f1[0] += ff * dpd1x[i];
        f1[1] += ff * dpd1y[i];
        f1[2] += ff * dpd1z[i];
        
        float* f2 = f + atoms[i][1] * 3;
        f2[0] += ff * dpd2x[i];
        f2[1] += ff * dpd2y[i];
        f2[2] += ff * dpd2z[i];
        
        float* f3 = f + atoms[i][2] * 3;
        f3[0] += ff * dpd3x[i];
        f3[1] += ff * dpd3y[i];
        f3[2] += ff * dpd3z[i];
        
        float* f4 = f + atoms[i][3] * 3;
        f4[0] += ff * dpd4x[i];
        f4[1] += ff * dpd4y[i];
        f4[2] += ff * dpd4z[i];
    }
}
