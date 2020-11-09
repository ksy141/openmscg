#include "pair_list.h"

#include <cmath>
#include <cstdio>

#undef NDEBUG
#include <cassert>

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
    int nx = static_cast<int>(ceil(pair->cut / pair->binsizex));
    int ny = static_cast<int>(ceil(pair->cut / pair->binsizey));
    int nz = static_cast<int>(ceil(pair->cut / pair->binsizez));
    
    
    
    int max_neigh = (nx*2+1) * (ny*2+1) * (nz*2+1);
    if(max_neigh<27) max_neigh = 27;
    
    neigh_bins = new int[max_neigh];
    sx = new int[max_neigh];
    sy = new int[max_neigh];
    sz = new int[max_neigh];
    
    int xme, yme, zme;
    pair->bin2offset(me, &xme, &yme, &zme);
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

PairList::PairList(float cut, float binsize, long maxpairs)
{
    this->cut = cut;
    this->cut_sq = cut * cut;
    this->binsize = binsize;
    this->maxpairs = maxpairs;
    
    maxbins = 0;
    binhead = links = ibins = 0;
    stencil = 0;
    
    ilist = new int[maxpairs];
    jlist = new int[maxpairs];
    tlist = new int[maxpairs];
    dxlist = new float[maxpairs];
    dylist = new float[maxpairs];
    dzlist = new float[maxpairs];
    drlist = new float[maxpairs];
}

PairList::~PairList()
{
    deallocate();
    if(stencil) delete [] stencil;
    
    if(binhead) delete [] binhead;
    if(ilist) delete [] ilist;
    if(jlist) delete [] jlist;
    if(tlist) delete [] tlist;
    if(dxlist) delete [] dxlist;
    if(dylist) delete [] dylist;
    if(dzlist) delete [] dzlist;
}

void PairList::allocate()
{
    links = new int [natoms];
    ibins = new int [natoms];
}

void PairList::deallocate()
{
    if(links) delete [] links;
    if(ibins) delete [] ibins;
}

void PairList::init(int* types, int natoms, int* exmap, int maxex)
{
    this->types = types;
    this->natoms = natoms;
    this->exmap = exmap;
    this->maxex = maxex;
    
    deallocate();
    allocate();
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

void PairList::setup_bins(vec3f box)
{
    for(int d=0; d<3; d++) 
    {
        this->box[d] = box[d];
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

void PairList::bin_atoms(vec3f *x)
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

inline bool not_excluded(int *ex, int target)
{
    int i=0; while(ex[i]!= -1)
    {
        if(ex[i]==target) return false; i += 1;
    } 
    return true;
}

void PairList::build(vec3f *x)
{
    bin_atoms(x);
    npairs = 0;
        
    for(int i=0; i<natoms; i++)
    {
        int *ex = 0;
        if (maxex) ex = exmap + maxex * i;
        
        double xi = x[i][0];
        double yi = x[i][1];
        double zi = x[i][2];
        
        int ibin = ibins[i];
        int j = binhead[ibin];
        
        while(j>=0)
        {
            if(j>i && (ex==0 || not_excluded(ex, j)))
            {
                double dx = x[j][0] - xi;
                double dy = x[j][1] - yi;
                double dz = x[j][2] - zi;
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
                if(ex==0 || not_excluded(ex, j))
                {
                    double dx = x[j][0] - xi + sx;
                    double dy = x[j][1] - yi + sy;
                    double dz = x[j][2] - zi + sz;
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
                
                j = links[j];
            }
        }
    }
    
    update_types(types);
}

void PairList::update_types(int *types)
{
    for(int i=0; i<npairs; i++)
        tlist[i] = pair_tid(types[ilist[i]], types[jlist[i]]);
}

