#include "table_pair_bspline.h"
#include "pair_list.h"

#include "std.h"

TablePairBSpline::TablePairBSpline(PairList *pair, int type_id, int order, double resolution, double xmin) : Table(), BSpline(order, resolution, xmin, pair->cut, 1)
{
    this->pair = pair;
    this->type_id = type_id;
    this->ncols = ncoeff;
}

TablePairBSpline::~TablePairBSpline()
{

}

void TablePairBSpline::compute()
{
    int natoms = pair->natoms;
    int* ilist = pair->ilist;
    int* jlist = pair->jlist;
    int* tlist = pair->tlist;
    float* dx = pair->dxlist;
    float* dy = pair->dylist;
    float* dz = pair->dzlist;
    float* dr = pair->drlist;

    for(int p=0; p<pair->npairs; p++) if(tlist[p]==type_id)
    {
        double *b;
        size_t istart;
        int nn;
        
        eval_coeffs(dr[p], &b, &istart, &nn);

        int atom_i = ilist[p];
        int atom_j = jlist[p];

        double *coeff_i, *coeff_j;

        coeff_i = coeff[atom_i*3];
        coeff_j = coeff[atom_j*3];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dx[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = coeff[atom_i*3+1];
        coeff_j = coeff[atom_j*3+1];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dy[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
 
        coeff_i = coeff[atom_i*3+2];
        coeff_j = coeff[atom_j*3+2];

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dz[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
    }
}

void TablePairBSpline::dump(double *output, double xmin, double dx, int n)
{
    eval(solution, output, xmin, dx, n);
}
