#include "model_angle_bspline.h"
#include "bond_list.h"
#include "defs.h"

ModelAngleBSpline::ModelAngleBSpline(double xmin, double xmax, double resolution, int order) : BSpline(order, resolution, xmin, xmax)
{
    nparam = ncoeff;
}

ModelAngleBSpline::~ModelAngleBSpline()
{

}

void ModelAngleBSpline::compute_fm()
{
    BondList *lst = (BondList*)(this->list);
    vec3i *atoms = lst->angle_atoms;
    int* types = lst->angle_types;
    
    float *dx1 = lst->dx1_angle;
    float *dy1 = lst->dy1_angle;
    float *dz1 = lst->dz1_angle;
    float *dx2 = lst->dx2_angle;
    float *dy2 = lst->dy2_angle;
    float *dz2 = lst->dz2_angle;
    float *a11 = lst->a11_angle;
    float *a12 = lst->a12_angle;
    float *a22 = lst->a22_angle;
    float *theta = lst->theta_angle;
    
    for(int i=0; i<lst->nangles; i++) if(types[i] == tid)
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

        int ia = atoms[i][0], ib = atoms[i][1], ic = atoms[i][2];
        double *coeff_i, *coeff_j, *coeff_k;

        coeff_i = dF + ia*3* nparam;
        coeff_j = dF + ib*3* nparam;
        coeff_k = dF + ic*3* nparam;
        
        for(int c=0; c<nn; c++)
        {
            double Bi = b[c];
            int pos = istart + c;
            coeff_i[pos] += Bi * f1x;
            coeff_j[pos] -= Bi * (f1x + f3x);
            coeff_k[pos] += Bi * f3x;
        }

        coeff_i = dF + (ia*3+1)* nparam;
        coeff_j = dF + (ib*3+1)* nparam;
        coeff_k = dF + (ic*3+1)* nparam;
        
        for(int c=0; c<nn; c++)
        {
            double Bi = b[c];
            int pos = istart + c;
            coeff_i[pos] += Bi * f1y;
            coeff_j[pos] -= Bi * (f1y + f3y);
            coeff_k[pos] += Bi * f3y;
        }

        coeff_i = dF + (ia*3+2)* nparam;
        coeff_j = dF + (ib*3+2)* nparam;
        coeff_k = dF + (ic*3+2)* nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c];
            int pos = istart + c;
            coeff_i[pos] += Bi * f1z;
            coeff_j[pos] -= Bi * (f1z + f3z);
            coeff_k[pos] += Bi * f3z;
        }
    }
}

void ModelAngleBSpline::compute_rem()
{
    for(int i=0; i<nparam; i++) dU[i] = 0.0;
    
    BondList *lst = (BondList*)(this->list);
    int* types = lst->angle_types;
    float *theta = lst->theta_angle;
    
    for(int i=0; i<lst->nangles; i++) if(types[i] == tid)
    {
        double *b;
        size_t istart;
        int nn;
        
        eval_coeffs(theta[i], &b, &istart, &nn);
        for(int c=0; c<nn; c++) dU[istart+c] += b[c];
    }
}

void ModelAngleBSpline::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}
