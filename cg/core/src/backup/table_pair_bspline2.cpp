#include "table_pair_bspline2.h"
#include "pair_list.h"

#include "std.h"

TablePairBSpline::TablePairBSpline(PairList *pair, int type_id, int order, double resolution, double xmin) : Table()
{
    this->type_id = type_id;
    this->resolution = resolution;
    this->xmin = xmin;
    this->pair = pair;
    this->xmax = pair->cut;
    int nbreak = static_cast<int>(ceil((xmax - xmin)/resolution)) + 1;

    this->order = order;
    bw = gsl_bspline_alloc(order, nbreak);
    gsl_bspline_knots_uniform(xmin, xmax, bw);

    ncols = order + nbreak - 2;
    B = gsl_vector_calloc(order);

    table_bvalue = 0;
    table_istart = 0;
    table_nn = 0;
    ddx = 0.0;
}

TablePairBSpline::~TablePairBSpline()
{
    gsl_bspline_free(bw);
    gsl_vector_free(B);

    if(table_bvalue)
    {
        delete [] table_bvalue;
        delete [] table_istart;
        delete [] table_nn;
    }
}

void TablePairBSpline::setup_cache(double dx_factor)
{
    ddx = dx_factor * resolution;
    int tsize = static_cast<int>(round((xmax - xmin) / ddx));
    table_bvalue = new double [tsize * order];
    table_istart = new int [tsize];
    table_nn = new int [tsize];

    for(size_t i=0; i<tsize; i++)
    {
        size_t istart, iend;
        double dx = xmin + ddx * i;
        gsl_bspline_eval_nonzero(dx, B, &istart, &iend, bw);
        int nn = iend - istart + 1;

        table_istart[i] = istart;
        table_nn[i] = nn;
        for(int j=0; j<nn; j++) table_bvalue[i*order+j] = B->data[j] / (-dx);
    }
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

    for(int p=0; p<pair->npairs; p++) if(tlist[p]==type_id && dr[p]>xmin && dr[p]<xmax)
    {
        size_t istart, iend;
        int nn;
        double *b;

        if(table_bvalue)
        {
            int it = static_cast<int>(round((dr[p] - xmin) / ddx));
            b = table_bvalue + it * order;
            istart = table_istart[it];
            nn = table_nn[it];
        }
        else
        {
            gsl_bspline_eval_nonzero(dr[p], B, &istart, &iend, bw);
            nn = iend - istart + 1;
            b = B->data;
            for(int c=0; c<=nn; c++) b[c] /= - dr[p];
        }

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
    //for(int i=0; i<ncols; i++) printf("%d %d %lf\n", i, ncols, solution[i]);
    
    gsl_vector *Bs = gsl_vector_calloc(ncols);
    gsl_vector *Cs = gsl_vector_calloc(ncols);

    for(int i=0; i<ncols; i++) gsl_vector_set(Cs, i, solution[i]);

    for(int i=0; i<n; i++)
    {
        gsl_bspline_eval(xmin + dx * i, Bs, bw);
        gsl_blas_ddot(Bs, Cs, output++);
    }

    gsl_vector_free(Bs);
    gsl_vector_free(Cs);
}
