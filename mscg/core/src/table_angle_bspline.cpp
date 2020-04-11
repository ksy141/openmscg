#include "table_angle_bspline.h"
#include "topology.h"
#include "bond_list.h"
#include "std.h"

TableAngleBSpline::TableAngleBSpline(BondList *lst, int type_id,
    int order, double resolution, double xmin, double xmax) : Table(), BSpline(order, resolution, xmin, xmax, 0)
{
    this->lst = lst;
    this->top = lst->top;
    this->type_id = type_id;
    this->ncols = ncoeff;
}

TableAngleBSpline::~TableAngleBSpline()
{

}

void TableAngleBSpline::compute()
{
    int *types = top->angl_types;
    int *atom1 = top->angl_atom1;
    int *atom2 = top->angl_atom2;
    int *atom3 = top->angl_atom3;
    
    float *dx1 = lst->dx1_angl;
    float *dy1 = lst->dy1_angl;
    float *dz1 = lst->dz1_angl;
    float *dx2 = lst->dx2_angl;
    float *dy2 = lst->dy2_angl;
    float *dz2 = lst->dz2_angl;
    float *a11 = lst->a11_angl;
    float *a12 = lst->a12_angl;
    float *a22 = lst->a22_angl;
    float *theta = lst->theta_angl;
    
    float r2d = 180 / 3.14159265359;
    
    for(int i=0; i<top->nangls; i++) if(types[i] == type_id)
    {         
        double *b;
        size_t istart;
        int nn;
        
        eval_coeffs(theta[i], &b, &istart, &nn);
        
        double f1x = a11[i]*dx1[i] + a12[i]*dx2[i];
        double f1y = a11[i]*dy1[i] + a12[i]*dy2[i];
        double f1z = a11[i]*dz1[i] + a12[i]*dz2[i];
        double f3x = a22[i]*dx2[i] + a12[i]*dx1[i];
        double f3y = a22[i]*dy2[i] + a12[i]*dy1[i];
        double f3z = a22[i]*dz2[i] + a12[i]*dz1[i];

        int ia = atom1[i], ib = atom2[i], ic = atom3[i];
        double *coeff_i, *coeff_j, *coeff_k;

        coeff_i = coeff[ia*3];
        coeff_j = coeff[ib*3];
        coeff_k = coeff[ic*3];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * r2d;
            int pos = istart + c;
            coeff_i[pos] += Bi * f1x;
            coeff_j[pos] -= Bi * (f1x + f3x);
            coeff_k[pos] += Bi * f3x;
        }

        coeff_i = coeff[ia*3+1];
        coeff_j = coeff[ib*3+1];
        coeff_k = coeff[ic*3+1];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * r2d;
            int pos = istart + c;
            coeff_i[pos] += Bi * f1y;
            coeff_j[pos] -= Bi * (f1y + f3y);
            coeff_k[pos] += Bi * f3y;
        }

        coeff_i = coeff[ia*3+2];
        coeff_j = coeff[ib*3+2];
        coeff_k = coeff[ic*3+2];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * r2d;
            int pos = istart + c;
            coeff_i[pos] += Bi * f1z;
            coeff_j[pos] -= Bi * (f1z + f3z);
            coeff_k[pos] += Bi * f3z;
        }
    }
}

void TableAngleBSpline::dump(double *output, double xmin, double dx, int n)
{
    eval(solution, output, xmin, dx, n);
}
