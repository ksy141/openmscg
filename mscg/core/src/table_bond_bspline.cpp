#include "table_bond_bspline.h"
#include "topology.h"
#include "bond_list.h"
#include "std.h"

TableBondBSpline::TableBondBSpline(BondList *lst, int type_id,
    int order, double resolution, double xmin, double xmax) : Table(), BSpline(order, resolution, xmin, xmax, 1)
{
    this->lst = lst;
    this->top = lst->top;
    this->type_id = type_id;
    this->ncols = ncoeff;
}

TableBondBSpline::~TableBondBSpline()
{
    
}

void TableBondBSpline::compute()
{
    int *types = top->bond_types;
    int *atom1 = top->bond_atom1;
    int *atom2 = top->bond_atom2;
    
    float *dx = lst->dx_bond;
    float *dy = lst->dy_bond;
    float *dz = lst->dz_bond;
    float *dr = lst->dr_bond;
        
    for(int i=0; i<top->nbonds; i++) if(types[i] == type_id)
    {         
        double *b;
        size_t istart;
        int nn;

        eval_coeffs(dr[i], &b, &istart, &nn);

        int ia = atom1[i], ib = atom2[i];
        double *coeff_i, *coeff_j;

        coeff_i = coeff[ia*3];
        coeff_j = coeff[ib*3];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dx[i];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = coeff[ia*3+1];
        coeff_j = coeff[ib*3+1];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dy[i];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = coeff[ia*3+2];
        coeff_j = coeff[ib*3+2];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dz[i];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
    }
}

void TableBondBSpline::dump(double *output, double xmin, double dx, int n)
{
    eval(solution, output, xmin, dx, n);
}
