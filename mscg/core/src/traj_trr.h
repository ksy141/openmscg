#ifndef TRAJ_TRR_H
#define TRAJ_TRR_H

#include "traj.h"
#include "xdrfile_trr.h"

class TrajTRR : public Traj
{
  public:
    
    t_trnheader header;
    XDRFILE *xd;
    
    TrajTRR(const char*, const char*);
    virtual ~TrajTRR();
    
    virtual int read_next_frame();
    virtual int write_frame();
    virtual void rewind();
};

#endif