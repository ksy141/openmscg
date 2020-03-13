#ifndef TRAJ_H
#define TRAJ_H

#include <cstdio>

typedef float Vec[3];

class Traj
{
  public:
    
    int status; /* 0 - ready, 1 - file error, 2 - format error, 3 - no data */
    int natoms;
    bool has_force;
    
    Vec box;
    Vec *x, *v, *f;
    
    Traj();
    virtual ~Traj();
    
    void pbc();
    virtual int read_next_frame() = 0;
    virtual void rewind() = 0;
};

#endif
