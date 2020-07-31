#include "traj.h"

Traj::Traj()
{
    natoms = maxatoms = 0;
    x = v = f = 0;
    type = 0;
    has_type = has_vel = has_force = false;
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
        maxatoms = natoms;
        type = new int[maxatoms];
        x = new Vec[maxatoms];
        v = new Vec[maxatoms];
        f = new Vec[maxatoms];
    }
}

void Traj::deallocate()
{
    if(type) delete [] type;
    if(x) delete [] x;
    if(v) delete [] v;
    if(f) delete [] f;
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
