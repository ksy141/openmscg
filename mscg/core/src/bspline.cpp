#include "bspline.h"

BSpline::BSpline(int order, double resolution, double xmin, double xmax, int scale_flag)
{
    //printf("BSpline %d, %lf, %lf, %lf\n", order, resolution, xmin, xmax);
    
    this->order = order;
    this->resolution = resolution;
    this->xmin = xmin;
    this->xmax = xmax;
    this->scale_flag = scale_flag;
    nbreak = static_cast<int>(ceil((xmax - xmin)/resolution)) + 1;

    bw = gsl_bspline_alloc(order, nbreak);
    gsl_bspline_knots_uniform(xmin, xmax, bw);

    ncoeff = order + nbreak - 2;
    B =  gsl_vector_calloc(order);
    
    table_bvalue = 0;
    table_istart = 0;
    table_nn = 0;
    ddx = 0.0;
    
    gsl_matrix *m = gsl_matrix_alloc(order, 2);
    gsl_matrix *mt = gsl_matrix_alloc(2, order);
    
    D0 = gsl_vector_calloc(order);
    D1 = gsl_vector_calloc(order);
    
    size_t iend;
    
    B0 = gsl_vector_calloc(order);
    D0 = gsl_vector_calloc(order);
    gsl_bspline_eval_nonzero(xmin, B0, &start0, &iend, bw);
    gsl_bspline_deriv_eval_nonzero(xmin, 1, m, &start0, &iend, bw);
    gsl_matrix_transpose_memcpy(mt, m);
    gsl_matrix_get_row(D0, mt, 1);
    
    B1 = gsl_vector_calloc(order);
    D1 = gsl_vector_calloc(order);
    gsl_bspline_eval_nonzero(xmax, B1, &start1, &iend, bw);
    gsl_bspline_deriv_eval_nonzero(xmax, 1, m, &start1, &iend, bw);
    gsl_matrix_transpose_memcpy(mt, m);
    gsl_matrix_get_row(D1, mt, 1);
        
    gsl_matrix_free(m);
    gsl_matrix_free(mt);
}

BSpline::~BSpline()
{
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_vector_free(D0);
    gsl_vector_free(B0);
    gsl_vector_free(D1);
    gsl_vector_free(B1);
    
    if(table_bvalue)
    {
        delete [] table_bvalue;
        delete [] table_istart;
        delete [] table_nn;
    }
}

void BSpline::setup_cache(double dx_factor)
{
    ddx = dx_factor * resolution;
    int tsize = static_cast<int>(round((xmax - xmin) / ddx)) + 1;
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
        
        if(scale_flag == 1) gsl_vector_scale(B, -1.0/dx);
        for(int j=0; j<nn; j++) table_bvalue[i*order+j] = B->data[j];
    }
}

void BSpline::eval_coeffs(double x, double **b, size_t *istart, int *nn)
{
    if(x<=xmin)
    {
        gsl_vector_memcpy(B, D0);
        gsl_vector_scale(B, (x-xmin));
        gsl_vector_add(B, B0);
        if(scale_flag) gsl_vector_scale(B, -1.0/x);
        (*istart) = start0;
        (*nn) = order;
        (*b) = B->data;
    }
    else if(x>=xmax)
    {
        gsl_vector_memcpy(B, D1);
        gsl_vector_scale(B, (x-xmax));
        gsl_vector_add(B, B1);
        if(scale_flag) gsl_vector_scale(B, -1.0/x);
        (*istart) = start1;
        (*nn) = order;
        (*b) = B->data;
    }
    else if(table_bvalue)
    {
        int it = static_cast<int>(round((x - xmin) / ddx));
        (*b) = table_bvalue + it * order;
        (*istart) = table_istart[it];
        (*nn) = table_nn[it];
    }
    else
    {
        size_t iend;
        gsl_bspline_eval_nonzero(x, B, istart, &iend, bw);        
        if(scale_flag==1) gsl_vector_scale(B, -1.0/x);
        
        (*nn) = iend - (*istart) + 1;
        (*b) = B->data;
    }
}

void BSpline::eval(double *input, double *output, double xmin, double dx, int n)
{
    gsl_vector *Bs = gsl_vector_calloc(ncoeff);
    gsl_vector *Cs = gsl_vector_calloc(ncoeff);
    
    for(int i=0; i<ncoeff; i++) gsl_vector_set(Cs, i, input[i]);

    for(int i=0; i<n; i++)
    {
        size_t istart;
        int nn;
        double *b;
        
        double x = xmin + dx * i;
        eval_coeffs(x, &b, &istart, &nn);
        
        gsl_vector_set_zero(Bs);
        for(int c=0; c<nn; c++) Bs->data[istart+c] = b[c];
        if(scale_flag) gsl_vector_scale(Bs, -x);
        gsl_blas_ddot(Bs, Cs, output++);
    }

    gsl_vector_free(Bs);
    gsl_vector_free(Cs);
}

void BSpline::eval(double *knots, double *in, double *out, int size)
{
    gsl_vector *Bs = gsl_vector_calloc(ncoeff);
    gsl_vector *Cs = gsl_vector_calloc(ncoeff);
    
    for(int i=0; i<ncoeff; i++) gsl_vector_set(Cs, i, knots[i]);

    for(int i=0; i<size; i++)
    {
        size_t istart;
        int nn;
        double *b;
        
        double x = in[i];
        eval_coeffs(x, &b, &istart, &nn);
        
        gsl_vector_set_zero(Bs);
        for(int c=0; c<nn; c++) Bs->data[istart+c] = b[c];
        if(scale_flag) gsl_vector_scale(Bs, -x);
        gsl_blas_ddot(Bs, Cs, out++);
    }

    gsl_vector_free(Bs);
    gsl_vector_free(Cs);
}


