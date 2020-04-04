#include "traj.h"

Traj::Traj()
{
    natoms = 0;
    x = v = f = 0;
}

Traj::~Traj()
{
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
        if(x[i][0]<0) x[i][0] += xprd;
        else if(x[i][0]>xprd) x[i][0] -= xprd;
        
        if(x[i][1]<0) x[i][1] += yprd;
        else if(x[i][1]>yprd) x[i][1] -= yprd;
        
        if(x[i][2]<0) x[i][2] += zprd;
        else if(x[i][2]>zprd) x[i][2] -= zprd;
    }
}