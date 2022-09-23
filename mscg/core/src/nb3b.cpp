#include "nb3b.h"

NB3B_SW::NB3B_SW(int i, int j, int k, float lambda, float cos0, float gamma_ij, float a_ij, float gamma_ik, float a_ik)
{
    this->type_i = i;
    this->type_j = j;
    this->type_k = k;
    this->lambda = lambda;
    this->cos0 = cos0;
    this->gamma_ij = gamma_ij;
    this->gamma_ik = gamma_ik;
}

NB3B_SW::~NB3B_SW()
{

}

void NB3B_SW::compute(PairList *plist, float *U, float *dU, float *f)
{
    int tid_ij = pair_tid(type_i, type_j);
    int tid_ik = pair_tid(type_i, type_k);

    int* tlist = plist->tlist;
    int* ilist = plist->ilist;
    int* jlist = plist->jlist;
    float* dx = plist->dxlist;
    float* dy = plist->dylist;
    float* dz = plist->dzlist;
    float* dr = plist->drlist;

    float dr1[3], dr2[3], fi[3], fj[3], fk[3];
    int n3 = 0;

    for(int i=0; i<plist->natoms; i++) if(plist->types[i] == type_i)
    {
        int j, k;
        int nneigh = plist->nneigh[i];
        int *neigh = plist->neigh_list[i];

        for(int n1=0; n1<nneigh; n1++)
        {
            int p1 = neigh[n1];

            if(tlist[p1] == tid_ij)
            {
                dr1[0] = dx[p1];
                dr1[1] = dy[p1];
                dr1[2] = dz[p1];
            }

            if(ilist[p1] == i) j = jlist[p1];
            else
            {
                j = ilist[p1];
                vector_reverse(dr1);
            }

            for(int n2=0; (n2<nneigh && n2!=n1); n2++)
            {
                int p2 = neigh[n2];

                if(tlist[p2] == tid_ik)
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

                    U[n3++] = compute_one(dr[p1], dr[p2], dr1, dr2, fi, fj, fk);

                    for(int d=0; d<3; d++)
                    {
                        f[i * 3 + d] += fi[d];
                        f[j * 3 + d] += fj[d];
                        f[k * 3 + d] += fk[d];
                    }
                }

            } // for(n2)
        } // for(n1)
    } // for(i)
}

float NB3B_SW::compute_one(float r1, float r2, float *dr1, float *dr2, float *fi, float *fj, float *fk)
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
    float dcs = cs - cos0;
    float dcssq = dcs * dcs;

    float dedxj = lambda * gamma_ij * rinv1a * rinv1a * rinv1 * e1e2 * dcssq;
    float dcdx = lambda * e1e2 * dcs * 2.0;
    float dcdxj1 = dcdx * cs * rinv1 * rinv1;
    float dcdxj2 = - dcdx * rinv_12;

    for(int d=0; d<3; d++)
    {
        fj[d] = (dedxj + dcdxj2) * dr1[d] + dcdxj1 * dr2[d];
        fk[d] = (dedxj + dcdxj2) * dr2[d] + dcdxj1 * dr1[d];
        fi[d] = - fj[d] - fk[d];
    }

    return lambda * dcssq * exp1 * exp2;
}
