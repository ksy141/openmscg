#include "matrix.h"
#include "table.h"
#include "std.h"

#include "mkl.h"
#define MEM_ALIGN 64

Matrix::Matrix()
{
    natoms  = 0;
    ncols   = 0;
    ntables = 0;
    
    matrix_coeff = 0;
    matrix_coeff_t = 0;
    matrix_cov = 0;
    vector_cov = 0;
    
    mkl_set_num_threads(1);
}

Matrix::~Matrix()
{
    if(matrix_coeff) mkl_free(matrix_coeff);
    if(matrix_coeff_t) mkl_free(matrix_coeff_t);
    if(matrix_cov) mkl_free(matrix_cov);
    if(vector_cov) mkl_free(vector_cov);
    if(vector_f) mkl_free(vector_f);
}

void Matrix::add_table(Table *tbl)
{
    tables[ntables++] = tbl;
    ncols += tbl->ncols;
}

void Matrix::setup(int natoms)
{
    this->natoms = natoms;
    size_t bytes;
    
    bytes = natoms * 3 * sizeof(double);
    vector_f = (double*)mkl_malloc(bytes, MEM_ALIGN);
    memset(vector_f, 0, bytes);
    
    bytes = ncols * ncols * sizeof(double);
    matrix_cov = (double*)mkl_malloc(bytes, MEM_ALIGN);
    memset(matrix_cov, 0, bytes);
    
    bytes = ncols * sizeof(double);
    vector_cov = (double*)mkl_malloc(bytes, MEM_ALIGN);
    memset(vector_cov, 0, bytes);
    
    bytes = natoms * 3 * ncols * sizeof(double);
    matrix_coeff = (double*)mkl_malloc(bytes, MEM_ALIGN);
    matrix_coeff_t = (double*)mkl_malloc(bytes, MEM_ALIGN);
    
    int icol = 0;
    
    for(int i=0; i<ntables; i++)
    {
        tables[i]->coeff = new double*[natoms * 3];
        
        for(int j=0; j<natoms*3; j++)
            tables[i]->coeff[j] = matrix_coeff + j * ncols + icol;
        
        tables[i]->solution = vector_cov + icol;
        icol += tables[i]->ncols;
    }
}

void Matrix::reset()
{
    size_t bytes = natoms * 3 * ncols * sizeof(double);
    memset(matrix_coeff, 0, bytes);
}

#define RANKK_OPERATOR 1

void Matrix::multiplyadd(float *F)
{    
    int k = natoms * 3;
    
    for(int i=0; i<k; i++) vector_f[i] = F[i];
        
    cblas_dgemv(CblasRowMajor, CblasTrans, natoms * 3, ncols,
        1.0, matrix_coeff, ncols, vector_f, 1, 1.0, vector_cov, 1);

#if RANKK_OPERATOR == 1
    double one = 1.0;
    mkl_domatcopy('R', 'T', k, ncols, 1.0, matrix_coeff, ncols, matrix_coeff_t, k);
    dsyrk("U", "T", &ncols, &k, &one, matrix_coeff_t, &k, &one, matrix_cov, &ncols);
    
#elif RANKK_OPERATOR == 2
    mkl_domatcopy('R', 'T', k, ncols, 1.0, matrix_coeff, ncols, matrix_coeff_t, k);
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, ncols, k,
        1.0, matrix_coeff_t, k, 1.0, matrix_cov, ncols);
    
#else
    size_t bytes = k * ncols * sizeof(double);
    memcpy(matrix_coeff_t, matrix_coeff, bytes);
    
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, ncols, ncols, k,
        1.0, matrix_coeff_t, ncols, matrix_coeff, ncols, 1.0, matrix_cov, ncols);

#endif
    
    /*
    int _t = 10;
    
    for(int i=0; i<_t; i++)
    {
        for(int j=0; j<_t; j++) printf(" %10.2lf", matrix_cov[i*ncols+j]);
        printf("\n");
    }
    printf("\n");
    print_vector(vector_cov, ncols);
    printf("\n");
    */
}

void Matrix::solve()
{
#if RANKK_OPERATOR > 0
    LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', ncols, 1, matrix_cov, ncols, vector_cov, ncols);
#else
    LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', ncols, ncols, 1, matrix_cov, ncols, vector_cov, 1);
#endif
    
    /*
    int _t = 10;
    
    for(int i=0; i<_t; i++)
    {
        for(int j=0; j<_t; j++) printf(" %10.2lf", matrix_cov[i*ncols+j]);
        printf("\n");
    }
    printf("\n");
    print_vector(vector_cov, ncols);
    printf("\n");
    */
}





