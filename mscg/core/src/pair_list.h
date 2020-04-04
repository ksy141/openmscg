#ifndef PAIR_LIST_H
#define PAIR_LIST_H

#include "traj.h"
#include "topology.h"

#define MAXPAIR 2000000

class Stencil
{
  public:
    
    Stencil();
    virtual ~Stencil();
    
    void setup(class PairList*, int);
    
    int n_neigh;
    int *neigh_bins, *sx, *sy, *sz;    
};

class PairList
{
  public:
    
    // settings
    
    Topology* top;
    
    int natoms;
    float cut, cut_sq;
    Vec box, hbox;
    
    // verlet list
    
    int nbinx, nbiny, nbinz, nbins, maxbins;
    float binsize, binsizex, binsizey, binsizez;
    float bininvx, bininvy, bininvz;
    
    int *binhead;
    int *links, *ibins;
    Stencil *stencil;
    
    // pair list
    
    long npairs;   
    int *ilist, *jlist, *tlist;
    float *dxlist, *dylist, *dzlist, *drlist;
    
    // functions
    
    PairList(Topology*, int);
    virtual ~PairList();
    
    int offset2bin(int, int, int);
    void bin2offset(int, int*, int*, int*);
    
    void init(float, float);
    void bin_atoms(Vec*);
    
    // api
    
    void setup_bins(Traj*);
    void build(Traj*, bool reset_bins = false);
    void build_brutal(Traj*);
};

#endif
