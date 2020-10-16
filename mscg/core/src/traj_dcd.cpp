#include "traj_dcd.h"

TrajDCD::TrajDCD(const char* filename, const char* mode) : Traj()
{
    if(mode[0]=='w' || mode[0]=='a')
    {
        status = 1;
        return;
    }
    
    nread = 0;
    dcd = open_dcd_read(filename, "dcd", &natoms, &nsets);
    allocate();
    timestep.coords = &(x[0][0]);
    status = 0;
}

TrajDCD::~TrajDCD()
{
    close_file_read(dcd);
}

void TrajDCD::rewind()
{
    nread = 0;
    dcd_rewind(dcd);
}

int TrajDCD::read_next_frame()
{
    if(status) return 1;
    if(nread>=dcd->nsets) return 1;
    if(read_next_timestep(dcd, natoms, &timestep)) return 1;
    
    box[0] = timestep.A;
    box[1] = timestep.B;
    box[2] = timestep.C;
    pbc();
    return 0;
}

int TrajDCD::write_frame()
{
    return 1;
}
