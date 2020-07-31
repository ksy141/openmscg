#include "pair_list.h"

#include <cmath>
#include <cstdio>

#undef NDEBUG
#include <cassert>

inline bool not_special(int *specials, int nspecial, int target)
{
    for(int i=0; i<nspecial; i++) if (specials[i] == target) return false;
    return true;
}

Stencil::Stencil()
{
    n_neigh = 0;
    neigh_bins = sx = sy = sz = 0;
}

Stencil::~Stencil()
{
    if(neigh_bins) delete [] neigh_bins;
    if(sx) delete [] sx;
    if(sy) delete [] sy;
    if(sz) delete [] sz;
}

void Stencil::setup(PairList *pair, int me)
{
    int max_neigh = pair->nbins<27?27:pair->nbins;
    
    neigh_bins = new int[max_neigh];
    sx = new int[max_neigh];
    sy = new int[max_neigh];
    sz = new int[max_neigh];
    
    int xme, yme, zme;
    pair->bin2offset(me, &xme, &yme, &zme);
    
    int nx = static_cast<int>(ceil(pair->cut / pair->binsizex));
    int ny = static_cast<int>(ceil(pair->cut / pair->binsizey));
    int nz = static_cast<int>(ceil(pair->cut / pair->binsizez));
    
    n_neigh = 0;
    
    for(int ix=xme; ix<=xme+nx; ix++)
        for(int iy=yme-ny; iy<=yme+ny; iy++)
            for(int iz=zme-nz; iz<=zme+nz; iz++)
                if(ix>xme || iy>yme || (iy==yme && iz>zme))
                {
                    assert((n_neigh<max_neigh) && "bin_setup() overflow: please check box-size, bin-size and cutoff.");
                    
                    if(ix>=pair->nbinx) sx[n_neigh] = 1;
                    else sx[n_neigh] = 0;
                    
                    if(iy>=pair->nbiny) sy[n_neigh] = 1;
                    else if(iy<0) sy[n_neigh] = -1;
                    else sy[n_neigh] = 0;
                    
                    if(iz>=pair->nbinz) sz[n_neigh] = 1;
                    else if(iz<0) sz[n_neigh] = -1;
                    else sz[n_neigh] = 0;
                    
                    int bin_neigh = pair->offset2bin(
                        ix - sx[n_neigh] * pair->nbinx,
                        iy - sy[n_neigh] * pair->nbiny,
                        iz - sz[n_neigh] * pair->nbinz);
                    
                    neigh_bins[n_neigh++] = bin_neigh;
                }
}

PairList::PairList(Topology *top, int natoms)
{
    maxbins = 0;
    binhead = links = ibins = 0;
    stencil = 0;
    
    this->top = top;
    this->natoms = natoms;
    
    links = new int [natoms];
    ibins = new int [natoms];
}

void PairList::init(float cut, float binsize)
{    
    this->cut = cut;
    cut_sq = cut * cut;
    
    ilist = new int[MAXPAIR];
    jlist = new int[MAXPAIR];
    tlist = new int[MAXPAIR];
    dxlist = new float[MAXPAIR];
    dylist = new float[MAXPAIR];
    dzlist = new float[MAXPAIR];
    drlist = new float[MAXPAIR];
    
    this->binsize = binsize;
}

PairList::~PairList()
{
    if(links) delete [] links;
    if(ibins) delete [] ibins;
    if(binhead) delete [] binhead;
    if(stencil) delete [] stencil;
    if(ilist) delete [] ilist;
    if(jlist) delete [] jlist;
    if(tlist) delete [] tlist;
    if(dxlist) delete [] dxlist;
    if(dylist) delete [] dylist;
    if(dzlist) delete [] dzlist;
}

int PairList::offset2bin(int x, int y, int z)
{
    return x * nbiny * nbinz + y * nbinz + z;
}

void PairList::bin2offset(int b, int *x, int *y, int *z)
{
    (*x) = b / (nbiny * nbinz);
    b %= nbiny * nbinz;
    
    (*y) = b / nbinz;
    (*z) = b % nbinz;
}

void PairList::setup_bins(Traj* traj)
{
    for(int d=0; d<3; d++) 
    {
        box[d] = traj->box[d];
        hbox[d] = 0.5 * box[d];
    }
    
    float binsizeinv = 1.0 / binsize;
    
    nbinx = static_cast<int>(box[0] * binsizeinv);
    nbiny = static_cast<int>(box[1] * binsizeinv);
    nbinz = static_cast<int>(box[2] * binsizeinv);
    
    binsizex = box[0]/nbinx;
    binsizey = box[1]/nbiny;
    binsizez = box[2]/nbinz;

    bininvx = 1.0 / binsizex;
    bininvy = 1.0 / binsizey;
    bininvz = 1.0 / binsizez;
    
    nbins = nbinx * nbiny * nbinz;
    //printf("bins: %d x %d x %d = %d\n", nbinx, nbiny, nbinz, nbins);
    
    if(nbins > maxbins)
    {
        maxbins = nbins;
        
        if(binhead) delete [] binhead;
        binhead = new int [maxbins];
                
        if(stencil) delete [] stencil;
        stencil = new Stencil[maxbins];
        
        for(int i=0; i<maxbins; i++)
            stencil[i].setup(this, i);
    }
}

void PairList::bin_atoms(Vec *x)
{
    for(int i=0; i<nbins; i++) binhead[i] = -1;
    
    for(int i=0; i<natoms; i++)
    {
        int bx = static_cast<int>(x[i][0] * bininvx);
        int by = static_cast<int>(x[i][1] * bininvy);
        int bz = static_cast<int>(x[i][2] * bininvz);
        
        int bin = offset2bin(bx, by, bz);
        if(bin<0 || bin>=nbins) continue;
        
        ibins[i] = bin;
        links[i] = binhead[bin];
        binhead[bin] = i;
    }
}

inline int PAIRTYPE(int i, int j, int n)
{
    return i<j?(i*n+j):(j*n+i);
}

void PairList::build(Traj* traj, bool reset_bins)
{
    if(reset_bins) setup_bins(traj);
    
    Vec *x = traj->x;
    bin_atoms(x);
    npairs = 0;
    
    int ntype = top->ntypes_atom;
    int *types = top->atom_types;
    int *nspecials = top->nspecials;
    int **special_pairs = top->special_pairs;
    
    for(int i=0; i<natoms; i++)
    {
        int *specials = special_pairs[i];
        int nspecial = nspecials[i];
        
        double xi = x[i][0];
        double yi = x[i][1];
        double zi = x[i][2];
        
        int ibin = ibins[i];
        int j = binhead[ibin];
        
        while(j>=0)
        {
            if(j>i && not_special(specials, nspecial, j))
            {
                double dx = x[j][0] - xi;
                double dy = x[j][1] - yi;
                double dz = x[j][2] - zi;
                double r2 = dx*dx + dy*dy + dz*dz;
                
                if(r2 < cut_sq) 
                {
                    ilist[npairs] = i;
                    jlist[npairs] = j;
                    tlist[npairs] = PAIRTYPE(types[i], types[j], ntype);
                    dxlist[npairs] = dx;
                    dylist[npairs] = dy;
                    dzlist[npairs] = dz;
                    drlist[npairs] = sqrt(r2);
                    npairs++;
                }
            }
            
            j = links[j];
        }
                
        for(int si=0; si<stencil[ibin].n_neigh; si++)
        {
            int jbin = stencil[ibin].neigh_bins[si];
            
            double sx = box[0] * stencil[ibin].sx[si];
            double sy = box[1] * stencil[ibin].sy[si];
            double sz = box[2] * stencil[ibin].sz[si];
            
            int j = binhead[jbin];
            
            while(j>=0)
            {
                if(not_special(specials, nspecial, j))
                {
                    double dx = x[j][0] - xi + sx;
                    double dy = x[j][1] - yi + sy;
                    double dz = x[j][2] - zi + sz;
                    double r2 = dx*dx + dy*dy + dz*dz;

                    if(r2 < cut_sq) 
                    {
                        ilist[npairs] = i;
                        jlist[npairs] = j;
                        tlist[npairs] = PAIRTYPE(types[i], types[j], ntype);
                        dxlist[npairs] = dx;
                        dylist[npairs] = dy;
                        dzlist[npairs] = dz;
                        drlist[npairs] = sqrt(r2);
                        npairs++;
                    }
                }
                
                j = links[j];
            }
        }
    }
    
    // printf("npairs: %ld\n", npairs);
    // for(int i=0; i<npairs; i++)
    //     printf("%6d %6d %6d %10lf\n", tlist[i], ilist[i], jlist[i], drlist[i]);
}

void PairList::update_types(int ntype, int * types)
{
    for(int i=0; i<npairs; i++)
        tlist[i] = PAIRTYPE(types[ilist[i]], types[jlist[i]], ntype);
}

void PairList::build_brutal(Traj* traj)
{
    Vec *x = traj->x;
    int npairs = 0;
    
    for(int i=0; i<natoms; i++)
    {
        double xi = x[i][0];
        double yi = x[i][1];
        double zi = x[i][2];
        
        for(int j=i+1; j<natoms; j++)
        {
            double dx = x[j][0] - xi;
            double dy = x[j][1] - yi;
            double dz = x[j][2] - zi;
            
            if(dx<-hbox[0]) dx+= box[0];
            else if(dx>hbox[0]) dx -= box[0];
            
            if(dy<-hbox[1]) dy+= box[1];
            else if(dy>hbox[1]) dy -= box[1];
            
            if(dz<-hbox[2]) dz+= box[2];
            else if(dz>hbox[2]) dz -= box[2];
            
            double r2 = dx*dx + dy*dy + dz*dz;

            if(r2 < cut_sq) 
            {
                ilist[npairs] = i;
                jlist[npairs] = j;
                dxlist[npairs] = dx;
                dylist[npairs] = dy;
                dzlist[npairs] = dz;
                drlist[npairs] = sqrt(r2);
                npairs++;
            }
        }
    }
}
