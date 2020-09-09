#include "traj_lammps.h"

TrajLAMMPS::TrajLAMMPS(const char* filename, const char* mode) : Traj()
{
    status = 1;
    fp = fopen((char*)filename, mode);
    if(!fp) return;
    
    if(mode[0]=='w' || mode[0]=='a')
    {
        status = 0;
        return;
    }
    
    status = 2;
    if(read_head()) return;
    rewind();
    
    parse_columns();
    if(cid==-1 || cx==-1 || cy==-1 || cz==-1) return;
    
    attrs['t'] = (ctype!=-1);
    attrs['q'] = (cq!=-1);
    attrs['x'] = (cx!=-1 && cy!=-1 && cz!=-1);
    attrs['v'] = (cvx!=-1 && cvy!=-1 && cvz!=-1);
    attrs['f'] = (cfx!=-1 && cfy!=-1 && cfz!=-1);
    
    allocate();
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
    for(int i=0; i<2; i++) GETLINE();
    sscanf(line, "%d", &(step));
    
    for(int i=0; i<2; i++) GETLINE();
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
    
    cx = cy = cz = cvx = cvy = cvz = cfx = cfy = cfz = cid = ctype = cq = -1;
    
    #define KEY(k,name) else if(strcmp(argv[i], #name)==0) c##k = i
    
    for(int i=0; i<argn; i++)
    {
        if(0);
        KEY(x,x);
        KEY(x,xu);
        KEY(y,y);
        KEY(y,yu);
        KEY(z,z);
        KEY(z,zu);
        KEY(vx,vx);
        KEY(vy,vy);
        KEY(vz,vz);
        KEY(fx,fx);
        KEY(fy,fy);
        KEY(fz,fz);
        KEY(id,id);
        KEY(type,type);
        KEY(q,q);
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
        
        if(attrs['t']) t[id] = atoi(argv[ctype]);
        if(attrs['q']) q[id] = atof(argv[cq]);
        
        if(attrs['v'])
        {
            v[id][0] = atof(argv[cvx]);
            v[id][1] = atof(argv[cvy]);
            v[id][2] = atof(argv[cvz]);
        }
        
        if(attrs['f'])
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

int TrajLAMMPS::write_frame()
{    
    fprintf(fp, "ITEM: TIMESTEP\n%d\n", step);
    fprintf(fp, "ITEM: NUMBER OF ATOMS\n%d\n", natoms);
    
    fprintf(fp, "ITEM: BOX BOUNDS pp pp pp\n");
    for(int dim=0; dim<3; dim++) fprintf(fp, "%e %e\n", 0.0, box[dim]);
    
    fprintf(fp, "ITEM: ATOMS id");
    if(attrs['t']) fprintf(fp, " type");
    if(attrs['q']) fprintf(fp, " q");
    fprintf(fp, " x y z");
    if(attrs['v']) fprintf(fp, " vx vy vz");
    if(attrs['f']) fprintf(fp, " fx fy fz");
    fprintf(fp, "\n");
    
    for(int i=0; i<natoms; i++)
    {
        fprintf(fp, "%d", i);
        if(attrs['t']) fprintf(fp, " %d", t[i]);
        if(attrs['q']) fprintf(fp, " %f", q[i]);
        fprintf(fp, " %f %f %f", x[i][0], x[i][1], x[i][2]);
        if(attrs['v']) fprintf(fp, " %f %f %f", v[i][0], v[i][1], v[i][2]);
        if(attrs['f']) fprintf(fp, " %f %f %f", f[i][0], f[i][1], f[i][2]);
        fprintf(fp, "\n");
    }
    
    return 0;
}






