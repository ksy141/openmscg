#include "model_nb3b_sw.h"
#include "pair_list.h"
#include "defs.h"
#include <cassert>

ModelNB3BSW::ModelNB3BSW(double gamma_ij, double a_ij, double gamma_ik, double a_ik, double theta0)
{
    this->gamma_ij = gamma_ij;
    this->a_ij = a_ij;

    this->gamma_ik = gamma_ik;
    this->a_ik = a_ik;

    this->cos0 = cos(theta0 / RAD2DEG);

    nparam = 1;
}

ModelNB3BSW::~ModelNB3BSW()
{
}

void ModelNB3BSW::setup_ex(int tid_ij, int tid_ik)
{
    this->tid_ij = tid_ij;
    this->tid_ik = tid_ik;
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

                        //printf("3BODY %d-%d-%d r1=%f r2=%f\n", j, i, k, dr[p1], dr[p2]);
                        compute_fm_one(i, j, k, dr1, dr2, dr[p1], dr[p2]);
                    }

                } // for(n2)
            }
        } // for(n1)
    } // for(i)
}

void ModelNB3BSW::compute_fm_one(int i, int j, int k, float *dr1, float *dr2, float r1, float r2)
{
    /* -----------------------------------------------------------------------
     * cos[t] = cos(theta) - cos(theta0)
     * exp[ij] = exp[gamma_ij/(r_ij-a_ij)]
     * exp[ik] = exp[gamma_ik/(r_ik-a_ik)]
     * U = lambda * cos[t]^2 * exp[ij] * exp[ik]
     *
     * d(cos[t]^2)/dx_j = 2.0 [cos(theta)-cos(theta0)] * [d(cos(theta))/dx_j]
     * d(cos(theta))/dx_j = (x_k-x_i)/(r_ij*r_ik) - (x_j-x_i)*cos(theta)/r_ij^2
     * d(exp[ij])/dx_j = - exp[ij] * gamma_ij * (x_j-x_i) / (r_ij-a_ij)^2 / r_ij
     *
     * dU/dx_j = lambda * { d(cos[t]^2)/dx_j * exp[ij] * exp[ik] +
     *                      cos[t]^2 * [d(exp[ij])/dx_j] * exp[ik] }
     *
     * dU/dx_j / lambda = {
     *   2.0 * cos[t] * exp[ij] * exp[ik] *
     *     [ (x_k-x_i)/(r_ij*r_ik) - (x_j-x_i)*cos(theta)/r_ij^2 ] +
     *   cos[t]^2 * exp[ij] * exp[ik] *
     *     [ - exp[ij] * gamma_ij  * (x_j-x_i) / (r_ij-a_ij)^2 / r_ij ] }
     *
     * F[x_i] = - dU / dx_j
     * --------------------------------------------------------------------- */

    float rinv1 = 1.0 / r1;
    float rinv1a = 1.0 / (r1 - a_ij);
    float exp1 = exp(gamma_ij * rinv1a);

    float rinv2 = 1.0 / r2;
    float rinv2a = 1.0 / (r2 - a_ik);
    float exp2 = exp(gamma_ik * rinv2a);

    float e1e2 = exp1 * exp2;
    float rinv_12 = 1.0 / (r1 * r2);
    float cs = vector_dot(dr1, dr2) * rinv_12;

    float dedxj = gamma_ij * rinv1a * rinv1a * rinv1;
    float dedxk = gamma_ik * rinv2a * rinv2a * rinv2;
    float dcdxj = cs * rinv1 * rinv1;
    float dcdxk = cs * rinv2 * rinv2;
    float dcdx2 = - rinv_12;

    for(int d=0; d<3; d++)
    {
        double dcs = cs - cos0;
        double Bi = dcs * dcs * e1e2;
        float fj = Bi * dedxj * dr1[d];
        float fk = Bi * dedxk * dr2[d];

        Bi =  2.0 * (cs - cos0) * e1e2;
        fj += Bi * (dcdxj * dr1[d] + dcdx2 * dr2[d]);
        fk += Bi * (dcdxk * dr2[d] + dcdx2 * dr1[d]);

        dF[(j * 3 + d) * nparam] += fj;
        dF[(k * 3 + d) * nparam] += fk;
        dF[(i * 3 + d) * nparam] -= (fj + fk);
    }
}

void ModelNB3BSW::compute_rem()
{
    assert(0 && "Not implemented yet!");
}

void ModelNB3BSW::get_table(double *params, double* in, double* out, int size)
{
    out[0] = params[0];
}
