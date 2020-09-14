#ifndef MODEL_H
#define MODEL_H

class Model
{
  public:    
    int nparam;
    int tid;
    void *list;

    double *dU;
    double *dF;
    
    Model(int, void*, double*, double*);
    virtual ~Model();
    
    virtual void compute_fm() = 0;
    virtual void compute_rem() = 0;
    virtual void get_table(double *params, double* in, double* out, int size) = 0;
};

#endif
