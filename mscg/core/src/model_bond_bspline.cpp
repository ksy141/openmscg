#include "model_bond_bspline.h"
#include "bond_list.h"
#include "defs.h"

ModelBondBSpline::ModelBondBSpline(double xmin, double xmax, double resolution, int order, int tid, void *list, double *dF, double *dU) : BSpline(order, resolution, xmin, xmax), Model(tid, list, dF, dU)
{
    nparam = ncoeff;
}

ModelBondBSpline::~ModelBondBSpline()
{
}

void ModelBondBSpline::compute_fm()
{
    BondList *lst = (BondList*)(this->list);
    vec2i *atoms = lst->bond_atoms;
    int* types = lst->bond_types;
    float *dx = lst->dx_bond;
    float *dy = lst->dy_bond;
    float *dz = lst->dz_bond;
    float *dr = lst->dr_bond;
        
    for(int i=0; i<lst->nbonds; i++) if(types[i] == tid)
    {         
        double *b;
        size_t istart;
        int nn;

        eval_coeffs(dr[i], &b, &istart, &nn);
        double inv = -1.0 / dr[i];
        
        int ia = atoms[i][0], ib = atoms[i][1];
        double *coeff_i, *coeff_j;

        coeff_i = dF + ia*3* nparam;
        coeff_j = dF + ib*3* nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dx[i] * inv;
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = dF + (ia*3+1)* nparam;
        coeff_j = dF + (ib*3+1)* nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dy[i] * inv;
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = dF + (ia*3+2)* nparam;
        coeff_j = dF + (ib*3+2)* nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dz[i] * inv;
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
    }
}

void ModelBondBSpline::compute_rem()
{
    
}

void ModelBondBSpline::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}