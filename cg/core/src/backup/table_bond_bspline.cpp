#include "table_bond_bspline.h"
#include "topology.h"
#include "bond_list.h"
#include "std.h"

TableBondBSpline::TableBondBSpline(BondList *lst, int type_id,
    int order, double resolution, double xmin, double xmax) : Table()
{
    this->lst = lst;
    this->top = lst->top;
    this->type_id = type_id;
    this->order = order;
    this->resolution = resolution;
    this->xmin = xmin;
    this->xmax = xmax;

    int nbreak = static_cast<int>(ceil((xmax - xmin)/resolution)) + 1;
    bw = gsl_bspline_alloc(order, nbreak);
    gsl_bspline_knots_uniform(xmin, xmax, bw);

    ncols = order + nbreak - 2;
    B = gsl_vector_calloc(order);

    table_bvalue = 0;
    table_istart = 0;
    table_nn = 0;
    ddx = 0.0;
}

TableBondBSpline::~TableBondBSpline()
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

void TableBondBSpline::setup_cache(double dx_factor)
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

void TableBondBSpline::compute()
{
    int *types = top->bond_types;
    int *atom1 = top->bond_atom1;
    int *atom2 = top->bond_atom2;
    
    float *dx = lst->dx_bond;
    float *dy = lst->dy_bond;
    float *dz = lst->dz_bond;
    float *dr = lst->dr_bond;
    
    for(int i=0; i<top->nbonds; i++)
    {
        if(types[i] == type_id)
        {         
            size_t istart, iend;
            int nn;
            double *b;
            
            
            if(dr[i]<xmin || dr[i]>xmax)
            {
                //printf("Error bond range: %f out of (%f, %f)\n", dr[i], xmin, xmax);
                continue;
            }
            
            
            if(table_bvalue)
            {
                int it = static_cast<int>(round((dr[i] - xmin) / ddx));
                b = table_bvalue + it * order;
                istart = table_istart[it];
                nn = table_nn[it];
            }
            else
            {
                gsl_bspline_eval_nonzero(dr[i], B, &istart, &iend, bw);
                nn = iend - istart + 1;
                b = B->data;
                for(int c=0; c<=nn; c++) b[c] /= - dr[i];
            }
            
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
}

void TableBondBSpline::dump(double *output, double xmin, double dx, int n)
{
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
