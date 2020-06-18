#ifndef DELTA_H
#define DELTA_H

class Delta
{
  public:
    
    double **coeff;
    double *solution;
    
    int ncols;
    
    Delta();
    virtual ~Delta();
    
    virtual void compute() = 0;
};

#endif
