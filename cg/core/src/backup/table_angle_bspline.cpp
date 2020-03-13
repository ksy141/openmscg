#include "table_angle_bspline.h"
#include "topology.h"
#include "bond_list.h"
#include "std.h"

TableAngleBSpline::TableAngleBSpline(BondList *lst, int type_id,
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

TableAngleBSpline::~TableAngleBSpline()
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

void TableAngleBSpline::setup_cache(double dx_factor)
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

void TableAngleBSpline::compute()
{
    int *types = top->angl_types;
    int *atom1 = top->angl_atom1;
    int *atom2 = top->angl_atom2;
    int *atom3 = top->angl_atom3;
    
    float *dx1 = lst->dx1_angl;
    float *dy1 = lst->dy1_angl;
    float *dz1 = lst->dz1_angl;
    float *dx2 = lst->dx2_angl;
    float *dy2 = lst->dy2_angl;
    float *dz2 = lst->dz2_angl;
    float *a11 = lst->a11_angl;
    float *a12 = lst->a12_angl;
    float *a22 = lst->a22_angl;
    float *theta = lst->theta_angl;
    
    for(int i=0; i<top->nangls; i++)
    {
        if(types[i] == type_id)
        {         
            size_t istart, iend;
            int nn;
            double *b;
            
            if(theta[i]<xmin || theta[i]>xmax)
            {
                //printf("Error bond range: %f out of (%f, %f)\n", theta[i], xmin, xmax);
                continue;
            }
            
            if(table_bvalue)
            {
                int it = static_cast<int>(floor((theta[i] - xmin) / ddx));
                b = table_bvalue + it * order;
                istart = table_istart[it];
                nn = table_nn[it];
            }
            else
            {
                gsl_bspline_eval_nonzero(theta[i], B, &istart, &iend, bw);
                nn = iend - istart + 1;
                b = B->data;
                for(int c=0; c<=nn; c++) b[c] /= - 1.0;
            }
            
            double f1x = a11[i]*dx1[i] + a12[i]*dx2[i];
            double f1y = a11[i]*dy1[i] + a12[i]*dy2[i];
            double f1z = a11[i]*dz1[i] + a12[i]*dz2[i];
            double f3x = a22[i]*dx2[i] + a12[i]*dx1[i];
            double f3y = a22[i]*dy2[i] + a12[i]*dy1[i];
            double f3z = a22[i]*dz2[i] + a12[i]*dz1[i];
            
            int ia = atom1[i], ib = atom2[i], ic = atom3[i];
            double *coeff_i, *coeff_j, *coeff_k;
            
            coeff_i = coeff[ia*3];
            coeff_j = coeff[ib*3];
            coeff_k = coeff[ic*3];
            
            for(int c=0; c<nn; c++)
            {
                double Bi = b[c];
                int pos = istart + c;
                coeff_i[pos] += Bi * f1x;
                coeff_j[pos] -= Bi * (f1x - f3x);
                coeff_k[pos] += Bi * f3x;
            }

            coeff_i = coeff[ia*3+1];
            coeff_j = coeff[ib*3+1];
            coeff_k = coeff[ic*3+1];

            for(int c=0; c<nn; c++)
            {
                double Bi = b[c];
                int pos = istart + c;
                coeff_i[pos] += Bi * f1y;
                coeff_j[pos] -= Bi * (f1y - f3y);
                coeff_k[pos] += Bi * f3y;
            }

            coeff_i = coeff[ia*3+2];
            coeff_j = coeff[ib*3+2];
            coeff_k = coeff[ic*3+2];

            for(int c=0; c<nn; c++)
            {
                double Bi = b[c];
                int pos = istart + c;
                coeff_i[pos] += Bi * f1z;
                coeff_j[pos] -= Bi * (f1z - f3z);
                coeff_k[pos] += Bi * f3z;
            }
        }
    }
}

void TableAngleBSpline::dump(double *output, double xmin, double dx, int n)
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
