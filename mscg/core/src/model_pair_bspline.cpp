#include "model_pair_bspline.h"
#include "pair_list.h"
#include "defs.h"

ModelPairBSpline::ModelPairBSpline(double xmin, double xmax, double resolution, int order, int tid, void *list, double *dF, double *dU) : BSpline(order, resolution, xmin, xmax, 1), Model(tid, list, dF, dU)
{
    nparam = ncoeff;
}

ModelPairBSpline::~ModelPairBSpline()
{
}

void ModelPairBSpline::compute_fm()
{
    PairList *pair = (PairList*)(this->list);
    int* ilist = pair->ilist;
    int* jlist = pair->jlist;
    int* tlist = pair->tlist;
    float* dx = pair->dxlist;
    float* dy = pair->dylist;
    float* dz = pair->dzlist;
    float* dr = pair->drlist;

    for(int p=0; p<pair->npairs; p++) if(tlist[p]==tid)
    {
        double *b;
        size_t istart;
        int nn;
        
        eval_coeffs(dr[p], &b, &istart, &nn);

        int atom_i = ilist[p];
        int atom_j = jlist[p];

        double *coeff_i, *coeff_j;

        coeff_i = dF + atom_i*3 * nparam;
        coeff_j = dF + atom_j*3 * nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dx[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }

        coeff_i = dF + (atom_i*3+1) * nparam;
        coeff_j = dF + (atom_j*3+1) * nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dy[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
 
        coeff_i = dF + (atom_i*3+2) * nparam;
        coeff_j = dF + (atom_j*3+2) * nparam;

        for(int c=0; c<nn; c++)
        {
            double Bi = b[c] * dz[p];
            int pos = istart + c;
            coeff_i[pos] += Bi;
            coeff_j[pos] -= Bi;
        }
    }
}

void ModelPairBSpline::compute_rem()
{
    for(int i=0; i<nparam; i++) dU[i] = 0.0;
    
    PairList *pair = (PairList*)list;
    int* tlist = pair->tlist;
    float* dr = pair->drlist;
        
    for(int p=0; p<pair->npairs; p++) if(tlist[p]==tid)
    {
        double *b;
        size_t istart;
        int nn;
        
        eval_coeffs(dr[p], &b, &istart, &nn);
        for(int c=0; c<nn; c++) dU[istart+c] = b[c];
    }
}

void ModelPairBSpline::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}