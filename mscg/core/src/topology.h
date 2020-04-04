#ifndef TOPOLOGY_H
#define TOPOLOGY_H

class Topology
{
  public:
    int ntypes_atom;
    
    int natoms, nbonds, nangls, ndihes;
    int maxatoms, maxbonds, maxangls, maxdihes;
    
    int *atom_types;
    int *bond_types, *bond_atom1, *bond_atom2;
    
    int *angl_types, *angl_atom1, *angl_atom2, *angl_atom3;
    int *dihe_types, *dihe_atom1, *dihe_atom2, *dihe_atom3, *dihe_atom4;
    
    int *nspecials, **special_pairs;
    
    Topology(int);
    virtual ~Topology();
    
    void add_atoms(int, int*);
    void add_bonds(int, int*, int*, int*);
    void add_angls(int, int*, int*, int*, int*);
    void add_dihes(int, int*, int*, int*, int*, int*);
    void build_special(bool, bool, bool);
};

#endif