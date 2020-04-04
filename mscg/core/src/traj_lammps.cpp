#include "traj_lammps.h"

TrajLAMMPS::TrajLAMMPS(const char* filename) : Traj()
{
    status = 1;
    fp = fopen((char*)filename, "r");
    if(!fp) return;
    
    status = 2;
    if(read_head()) return;
    rewind();
    
    parse_columns();
    if(cid==-1 || cx==-1 || cy==-1 || cz==-1) return;
    if(cfx!=-1 && cfy!=-1 && cfz!=-1) has_force = true;
    else has_force = false;
    
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

TrajLAMMPS::~TrajLAMMPS()
{
    if(fp) fclose(fp);
}

void TrajLAMMPS::rewind()
{
    if(fp) fseek(fp, 0, SEEK_SET);
}

int TrajLAMMPS::read_next_frame()
{
    if(read_head()) return 1;
    if(read_body()) return 1;
    return 0;
}

#define GETLINE() {if(fgets(line,1000,fp)==0) return 1;}

int TrajLAMMPS::read_head()
{
    for(int i=0; i<4; i++) GETLINE();
    sscanf(line, "%d", &(natoms));
    
    GETLINE();
    
    for(int i=0; i<3; i++)
    {
        GETLINE();
        sscanf(line, "%e %e", &(boxlo[i]), &(box[i]));
        box[i] -= boxlo[i];
    }
    
    GETLINE();
    strcpy(columns, line);
    
    return 0;
}

int TrajLAMMPS::parse_columns()
{
    int argn = 0;
    char *argv[200];
    
    strtok(columns, " \r\n"); // ITEMS:
    strtok(NULL, " \r\n"); // ATOMS
    char *ch = strtok(NULL, " \r\n");
    
    while(ch)
    {
        argv[argn++] = ch;
        ch = strtok(NULL, " \r\n");
    }
    
    cx = cy = cz = cfx = cfy = cfz = cid = ctype = -1;
    
    #define KEY(k) else if(strcmp(argv[i], #k)==0) c##k = i
    
    for(int i=0; i<argn; i++)
    {
        if(0);
        KEY(x);
        KEY(y);
        KEY(z);
        KEY(fx);
        KEY(fy);
        KEY(fz);
        KEY(id);
        KEY(type);
    }
    
    return 0;
}

int TrajLAMMPS::read_body()
{
    for(int i=0; i<natoms; i++)
    {
        GETLINE();
        
        int argn = 0;
        char *argv[200];
        char *ch = strtok(line, " \r\n");

        while(ch)
        {
            argv[argn++] = ch;
            ch = strtok(NULL, " \r\n");
        }
        
        int id = atoi(argv[cid]) - 1;
        
        x[id][0] = atof(argv[cx]);
        x[id][1] = atof(argv[cy]);
        x[id][2] = atof(argv[cz]);
        
        if(has_force)
        {
            f[id][0] = atof(argv[cfx]);
            f[id][1] = atof(argv[cfy]);
            f[id][2] = atof(argv[cfz]);
        }
    }
    
    for(int i=0; i<natoms; i++)
    {
        x[i][0] -= boxlo[0];
        x[i][1] -= boxlo[1];
        x[i][2] -= boxlo[2];
    }
    
    pbc();
    
    return 0;
}








