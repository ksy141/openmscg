#include "traj.h"

Traj::Traj()
{
    natoms = maxatoms = 0;
    
    t = 0;
    q = 0;
    x = v = f = 0;
    
    attrs['t'] = false;
    attrs['q'] = false;
    attrs['x'] = true;
    attrs['v'] = false;
    attrs['f'] = false;
}

Traj::~Traj()
{
    deallocate();
}

void Traj::allocate()
{
    if(natoms>maxatoms) 
    {
        deallocate();
        maxatoms = natoms + 10000;
        t = new int[maxatoms];
        q = new float[maxatoms];
        x = new Vec[maxatoms];
        v = new Vec[maxatoms];
        f = new Vec[maxatoms];
    }
}

void Traj::deallocate()
{
    if(t) delete [] t;
    if(q) delete [] q;
    if(x) delete [] x;
    if(v) delete [] v;
    if(f) delete [] f;
    
    t = 0;
    q = 0;
    x = v = f = 0;
}

void Traj::pbc()
{
    float xprd = box[0];
    float yprd = box[1];
    float zprd = box[2];
    
    for(int i=0; i<natoms; i++)
    {
        while(x[i][0]<0) x[i][0] += xprd;
        while(x[i][0]>xprd) x[i][0] -= xprd;
        
        while(x[i][1]<0) x[i][1] += yprd;
        while(x[i][1]>yprd) x[i][1] -= yprd;
        
        while(x[i][2]<0) x[i][2] += zprd;
        while(x[i][2]>zprd) x[i][2] -= zprd;
    }
}
