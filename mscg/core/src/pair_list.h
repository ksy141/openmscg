#ifndef PAIR_LIST_H
#define PAIR_LIST_H

#include "defs.h"

inline int pair_tid(int i, int j)
{
    if(i>j) return i * (i + 1) / 2 + j;
    else return j * (j + 1) / 2 + i;
}

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

    long maxpairs;
    int maxtypeid;
    int* types;
    int natoms;
    int* exmap;
    int maxex;
    float cut, cut_sq;
    vec3f box, hbox;

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

    // neighbor list of atoms
    int maxneigh;
    int *nneigh, *neigh_page;
    int **neigh_list;

    // functions

    PairList(float cut, float binsize, long maxpairs = 2000000);
    virtual ~PairList();

    void allocate();
    void deallocate();

    int offset2bin(int, int, int);
    void bin2offset(int, int*, int*, int*);
    void bin_atoms(vec3f* x);

    // api

    void init(int *types, int natoms, int* exmap, int maxex);
    void setup_bins(vec3f box);
    void build(vec3f *x);
    void build_neighbors();
    void update_types(int*);
    int count_3b();
};

#endif
