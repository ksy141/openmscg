#ifndef MODEL_H
#define MODEL_H

class Model
{
  public:    
    int nparam;
    int type_id;
    void *list;
    double *dudl;
    
    Model(void*, int);
    virtual ~Model();
    
    virtual void compute_dfdl() = 0;
    virtual void compute_dudl() = 0;
    virtual void compute_etable(double *params, double* in, double* out, int size) = 0;
    virtual void compute_ftable(double *params, double* in, double* out, int size) = 0;
};

#endif
