#include "model_nb3b_sw_bspline.h"
#include "pair_list.h"
#include "defs.h"
#include <cassert>

ModelNB3BSW::ModelNB3BSW(double xmin, double xmax, double resolution, int order) : BSpline(order, resolution, xmin, xmax)
{
    nparam = ncoeff;
}

ModelNB3BSW::~ModelNB3BSW()
{
}

void ModelNB3BSW::setup_ex(int tid_ij, double gamma_ij, double a_ij,
    int tid_ik, double gamma_ik, double a_ik)
{
    this->tid_ij = tid_ij;
    this->gamma_ij = gamma_ij;
    this->a_ij = a_ij;

    this->tid_ik = tid_ik;
    this->gamma_ik = gamma_ik;
    this->a_ik = a_ik;
}

void ModelNB3BSW::compute_fm()
{
    PairList *plist = (PairList*)(this->list);
    int* tlist = plist->tlist;
    int* ilist = plist->ilist;
    int* jlist = plist->jlist;
    float* dx = plist->dxlist;
    float* dy = plist->dylist;
    float* dz = plist->dzlist;
    float* dr = plist->drlist;

    float dr1[3], dr2[3];

    for(int i=0; i<plist->natoms; i++) if(plist->types[i] == tid)
    {
        int j, k;
        int nneigh = plist->nneigh[i];
        int *neigh = plist->neigh_list[i];

        for(int n1=0; n1<nneigh; n1++)
        {
            int p1 = neigh[n1];

            if(tlist[p1] == tid_ij && dr[p1]<a_ij)
            {
                dr1[0] = dx[p1];
                dr1[1] = dy[p1];
                dr1[2] = dz[p1];

                if(ilist[p1] == i) j = jlist[p1];
                else
                {
                    j = ilist[p1];
                    vector_reverse(dr1);
                }

                for(int n2=0; n2<nneigh; n2++)
                {
                    if(n1 == n2) continue;
                    int p2 = neigh[n2];

                    if(tlist[p2] == tid_ik && dr[p2]<a_ik)
                    {
                        dr2[0] = dx[p2];
                        dr2[1] = dy[p2];
                        dr2[2] = dz[p2];

                        if(ilist[p2] == i) k = jlist[p2];
                        else
                        {
                            k = ilist[p2];
                            vector_reverse(dr2);
                        }

                        printf("3BODY %d-%d-%d r1=%f r2=%f\n", i, j, k, dr[p1], dr[p2]);
                        compute_fm_one(i, j, k, dr1, dr2, dr[p1], dr[p2]);
                    }

                } // for(n2)
            }
        } // for(n1)
    } // for(i)
}

void ModelNB3BSW::compute_fm_one(int i, int j, int k, float *dr1, float *dr2, float r1, float r2)
{
    float rinv1 = 1.0 / r1;
    float rinv1a = 1.0 / (r1 - a_ij);
    float exp1 = exp(gamma_ij * rinv1a);

    //float rinv2 = 1.0 / r2;
    float rinv2a = 1.0 / (r2 - a_ik);
    float exp2 = exp(gamma_ik * rinv2a);
    float e1e2 = exp1 * exp2;
    printf("EXP %f\n", e1e2);
    float rinv_12 = 1.0 / (r1 * r2);
    float cs = vector_dot(dr1, dr2) * rinv_12;

    size_t istart_coeff, istart_deriv;
    int nn_coeff, nn_deriv;
    double *b_coeff, *b_deriv;
    eval_coeffs(cs, &b_coeff, &istart_coeff, &nn_coeff);
    eval_derivs(cs, &b_deriv, &istart_deriv, &nn_deriv);

    float dedxj = gamma_ij * rinv1a * rinv1a * rinv1 * e1e2;
    float dcdxj1 = cs * rinv1 * rinv1 * e1e2;
    float dcdxj2 = - rinv_12 * e1e2;

    for(int d=0; d<3; d++)
    {
        double *coeff_i = dF + (i * 3 + d) * nparam;
        double *coeff_j = dF + (j * 3 + d) * nparam;
        double *coeff_k = dF + (k * 3 + d) * nparam;

        for(int c=0; c<nn_coeff; c++)
        {
            double Bi = b_coeff[c] * dedxj;
            int pos = istart_coeff + c;

            coeff_i[pos] -= Bi * (dr1[d] + dr2[d]);
            coeff_j[pos] += Bi * dr1[d];
            coeff_k[pos] += Bi * dr2[d];
        }

        for(int c=0; c<nn_deriv; c++)
        {
            double Bi = b_deriv[c];
            int pos = istart_deriv + c;

            coeff_i[pos] -= Bi * (dcdxj1 + dcdxj2) * (dr1[d] + dr2[d]);
            coeff_j[pos] += Bi * (dcdxj1 * dr1[d] + dcdxj2 * dr2[d]);
            coeff_k[pos] += Bi * (dcdxj1 * dr2[d] + dcdxj2 * dr1[d]);
        }
    }
}

void ModelNB3BSW::compute_rem()
{
    assert(0 && "Not implemented yet!");
}

void ModelNB3BSW::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}
