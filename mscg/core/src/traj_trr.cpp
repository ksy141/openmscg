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

TrajTRR::TrajTRR(const char* filename, const char* mode) : Traj()
{
    status = 1;
    xd = xdrfile_open((char*)filename, mode);
    if(!xd) return;
    
    if(mode[0]=='w' || mode[0]=='a')
    {
        status = 0;
        return;
    }
    
    status = 2;
    if(do_trnheader(xd, 1, &header)) return;
    rewind();
    
    natoms = header.natoms;
    attrs['v'] = (header.v_size > 0);
    attrs['f'] = (header.f_size > 0);
    
    allocate();
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
    matrix _box;
    float t, lambda;
    
    if(status) return 1;
    if(read_trr(xd, natoms, &step, &t, &lambda, _box, x, v, f)) return 1;
    
    for(int i=0; i<natoms; i++)
    {
        x[i][0] *= 10.0; x[i][1] *= 10.0; x[i][2] *= 10.0;
        if(attrs['v']) { v[i][0] *= 0.01; v[i][1] *= 0.01; v[i][2] *= 0.01; }
        if(attrs['f']) { f[i][0] /= 41.82; f[i][1] /= 41.82; f[i][2] /= 41.82; }
    }
    
    if(_box[0][0]>1.0E-6 && _box[1][1]>1.0E-6 && _box[2][2]>1.0E-6)
        for(int dim=0; dim<3; dim++) box[dim] = _box[dim][dim] * 10.0;
    
    if(box[0]>1.0E-6 && box[1]>1.0E-6 && box[2]>1.0E-6) pbc();
    return 0;
}

int TrajTRR::write_frame()
{
    matrix _box;
    float lambda = 1.0;
    float t = step;
    
    for(int dim=0; dim<3; dim++)
    {
        _box[dim][0] = _box[dim][1] = _box[dim][2] = 0.0;
        _box[dim][dim] = box[dim] * 0.1;
    }
    
    for(int i=0; i<natoms; i++)
    {
        x[i][0] *= 0.1; x[i][1] *= 0.1; x[i][2] *= 0.1;
        if(attrs['v']) { v[i][0] *= 100.0; v[i][1] *= 100.0; v[i][2] *= 100.0; }
        if(attrs['f']) { f[i][0] *= 41.82; f[i][1] *= 41.82; f[i][2] *= 41.82; }
    }
    
    return write_trr(xd, natoms, step, t, lambda, _box, x, attrs['v']?v:NULL, attrs['f']?f:NULL);
}
