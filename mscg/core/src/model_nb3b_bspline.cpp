#include "model_nb3b_bspline.h"
#include "pair_list.h"
#include "defs.h"
#include <cassert>

ModelNB3BBSpline::ModelNB3BBSpline(double xmin, double xmax, double resolution, int order) : BSpline(order, resolution, xmin, xmax)
{
    nparam = ncoeff;
}

ModelNB3BBSpline::~ModelNB3BBSpline()
{
}

void ModelNB3BBSpline::setup_ex(int tid_ij, double gamma_ij, double a_ij,
    int tid_ik, double gamma_ik, double a_ik)
{
    this->tid_ij = tid_ij;
    this->gamma_ij = gamma_ij;
    this->a_ij = a_ij;

    this->tid_ik = tid_ik;
    this->gamma_ik = gamma_ik;
    this->a_ik = a_ik;

    this->using_theta = 1;
}

void ModelNB3BBSpline::compute_fm()
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
                    if(n2 == n1)
                    {
                        if(tid_ij == tid_ik) break;
                        else continue;
                    }

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

                        //printf("3BODY %d-%d-%d r1=%f r2=%f\n", i, j, k, dr[p1], dr[p2]);
                        compute_fm_one(i, j, k, dr1, dr2, dr[p1], dr[p2]);
                    }

                } // for(n2)
            }
        } // for(n1)
    } // for(i)
}

void ModelNB3BBSpline::compute_fm_one(int i, int j, int k, float *dr1, float *dr2, float r1, float r2)
{
    float rinv1 = 1.0 / r1;
    float rinv1a = 1.0 / (r1 - a_ij);
    float exp1 = exp(gamma_ij * rinv1a);

    float rinv2 = 1.0 / r2;
    float rinv2a = 1.0 / (r2 - a_ik);
    float exp2 = exp(gamma_ik * rinv2a);

    float e1e2 = exp1 * exp2;
    float rinv_12 = 1.0 / (r1 * r2);
    float cs = vector_dot(dr1, dr2) * rinv_12;
    float spline_x = cs;

    float theta = 0.0;

    if(using_theta)
    {
        theta = acos(cs);
        spline_x = theta / M_PI * 180;

        if(spline_x<0.1 || spline_x>179.99999) printf("ABNORMAL %d-%d-%d=%d (%f %f %f) (%f %f %f)\n", i, j, k, theta, dr1[0], dr1[1], dr1[2], dr2[0], dr2[1], dr2[2]);
    }

    size_t istart_coeff, istart_deriv;
    int nn_coeff, nn_deriv;
    double *b_coeff, *b_deriv;


    eval_coeffs(spline_x, &b_coeff, &istart_coeff, &nn_coeff);
    eval_derivs(spline_x, &b_deriv, &istart_deriv, &nn_deriv);

    float dedxj = gamma_ij * rinv1a * rinv1a * rinv1 * e1e2;
    float dedxk = gamma_ik * rinv2a * rinv2a * rinv2 * e1e2;
    float dcdxj = cs * rinv1 * rinv1 * e1e2;
    float dcdxk = cs * rinv2 * rinv2 * e1e2;
    float dcdx2 = - rinv_12 * e1e2;

    for(int d=0; d<3; d++)
    {
        double *coeff_i = dF + (i * 3 + d) * nparam;
        double *coeff_j = dF + (j * 3 + d) * nparam;
        double *coeff_k = dF + (k * 3 + d) * nparam;


        float fj = dedxj * dr1[d];
        float fk = dedxk * dr2[d];
        float fi = - fj - fk;

        for(int c=0; c<nn_coeff; c++)
        {
            double Bi = b_coeff[c];
            int pos = istart_coeff + c;

            coeff_i[pos] += Bi * fi;
            coeff_j[pos] += Bi * fj;
            coeff_k[pos] += Bi * fk;
        }

        fj = dcdxj * dr1[d] + dcdx2 * dr2[d];
        fk = dcdxk * dr2[d] + dcdx2 * dr1[d];
        fi = - fj - fk;

        for(int c=0; c<nn_deriv; c++)
        {
            double Bi = b_deriv[c];
            if(using_theta) Bi *= -1.0 / sin(theta);

            int pos = istart_deriv + c;

            coeff_i[pos] += Bi * fi;
            coeff_j[pos] += Bi * fj;
            coeff_k[pos] += Bi * fk;
        }
    }
}

void ModelNB3BBSpline::compute_rem()
{
    assert(0 && "Not implemented yet!");
}

void ModelNB3BBSpline::get_table(double *params, double* in, double* out, int size)
{
    eval(params, in, out, size);
}
