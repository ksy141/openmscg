#ifndef TRAJ_H
#define TRAJ_H

#include <cstdio>

typedef float Vec[3];

class Traj
{
  public:
    
    int status; /* 0 - ready, 1 - file error, 2 - format error, 3 - no data */
    int natoms;
    int maxatoms;
    int step;
    bool has_type;
    bool has_vel;
    bool has_force;
    
    Vec box;
    Vec *x, *v, *f;
    int *type;
    
    Traj();
    virtual ~Traj();
    
    void allocate();
    void deallocate();
    void pbc();
    
    virtual void rewind() = 0;
    virtual int read_next_frame() = 0;
    virtual int write_frame() = 0;
};

#endif
