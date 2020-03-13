#ifndef MATRIX_H
#define MATRIX_H

#define MAXTABLE 1024

class Matrix
{
  public:
    
    Matrix();
    virtual ~Matrix();
    
    int natoms;
    int ncols;
    void setup(int);
    
    int ntables;
    class Table* tables[MAXTABLE];
    void add_table(class Table*);
    
    double *matrix_coeff, *matrix_coeff_t;
    double *matrix_cov;
    double *vector_cov, *vector_f;
    
    void reset();
    void multiplyadd(float*);
    void solve();
};

#endif
