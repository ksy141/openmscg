#ifndef TABLE_H
#define TABLE_H

class Table
{
  public:
    
    double **coeff;
    double *solution;
    
    int ncols;
    
    Table();
    virtual ~Table();
    
    virtual void compute() = 0;
    virtual void dump(double*, double, double, int) = 0;
};

#endif
