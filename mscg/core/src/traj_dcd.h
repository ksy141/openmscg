#ifndef TRAJ_DCD_H
#define TRAJ_DCD_H

#include "traj.h"
#include "dcd/dcdplugin.h"

class TrajDCD : public Traj
{
  public:
    dcdhandle *dcd;
    int nread, nsets;
    molfile_timestep_t timestep;
    
    TrajDCD(const char*, const char*);
    virtual ~TrajDCD();
    
    virtual int read_next_frame();
    virtual int write_frame();
    virtual void rewind();
};

#endif