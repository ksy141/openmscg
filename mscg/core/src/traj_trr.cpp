#include "traj_trr.h"

struct XDRFILE
{
    FILE *   fp;
    struct XDR *    xdr;
    char     mode;
    int *    buf1; 
    int      buf1size;
    int *    buf2; 
    int      buf2size; 
};

TrajTRR::TrajTRR(const char* filename) : Traj()
{
    status = 1;
    xd = xdrfile_open((char*)filename, "r");
    if(!xd) return;
    
    status = 2;
    if(do_trnheader(xd, 1, &header)) return;
    rewind();
    
    natoms = header.natoms;
    has_force = (header.f_size > 0);
    
    x = new Vec[natoms];
    v = new Vec[natoms];
    f = new Vec[natoms];
    
    status = 0;
    
    if(read_next_frame()) 
    {
        status = 3;
        return;
    }
    
    rewind();
}

TrajTRR::~TrajTRR()
{
    if(xd) xdrfile_close(xd);
}

void TrajTRR::rewind()
{
    if(xd) fseek(xd->fp, 0, SEEK_SET);
}

int TrajTRR::read_next_frame()
{
    int step;
    matrix _box;
    float t, lambda;
    
    if(status) return 1;
    if(read_trr(xd, natoms, &step, &t, &lambda, _box, x, v, f)) return 1;
    for(int dim=0; dim<3; dim++) box[dim] = _box[dim][dim];
    pbc();
    return 0;
}
