#ifndef TABLES_H
#define TABLES_H

class Tables {
  public:
    
    struct Table {
        float min, inc;
        int n;
        float *efac, *ffac;
    } *tbls;
    
    Tables(int);
    ~Tables();
    
    void compute(long n, int *tid, float *x, float *e, float *f);
};

#endif