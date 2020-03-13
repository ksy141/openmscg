#include "topology.h"
#include <cstdlib>
#include <cstdio>

Topology::Topology(int ntypes_atom)
{
    this->ntypes_atom = ntypes_atom;
    
    natoms = nbonds = nangls = ndihes = 0;
    maxatoms = maxbonds = maxangls = maxdihes = 0;
    atom_types = bond_types = angl_types = dihe_types = 0;
    
    bond_atom1 = bond_atom2 = 0;
    angl_atom1 = angl_atom2 = angl_atom3 = 0;
    dihe_atom1 = dihe_atom2 = dihe_atom3 = dihe_atom4 = 0;
    
    nspecials = 0;
    special_pairs = 0;
}

Topology::~Topology()
{
    if(natoms) { free(atom_types); }
    if(nbonds) { free(bond_types); free(bond_atom1); free(bond_atom2); }
    
    if(nangls) { 
        free(angl_types);
        free(angl_atom1);
        free(angl_atom2);
        free(angl_atom3);
    }
    
    if(ndihes) {
        free(dihe_types);
        free(dihe_atom1);
        free(dihe_atom2);
        free(dihe_atom3);
        free(dihe_atom4);
    }
    
    if(nspecials) 
    {
        delete [] nspecials;
        delete [] special_pairs[0];
        delete [] special_pairs;
    }
}

void Topology::add_atoms(int n, int *types)
{
    for(int i=0; i<n; i++)
    {
        if(natoms >= maxatoms) {
            maxatoms += 1000;
            atom_types = (int*)realloc(atom_types, sizeof(int) * maxatoms);
        }
        
        atom_types[natoms++] = types[i];
    }
}

void Topology::add_bonds(int n, int *types, int *atom1, int *atom2)
{
    for(int i=0; i<n; i++)
    {
        if(nbonds >= maxbonds)
        {
            maxbonds += 1000;
            bond_types = (int*)realloc(bond_types, sizeof(int) * maxbonds);
            bond_atom1 = (int*)realloc(bond_atom1, sizeof(int) * maxbonds);
            bond_atom2 = (int*)realloc(bond_atom2, sizeof(int) * maxbonds);
        }
        
        bond_types[nbonds] = types[i];
        bond_atom1[nbonds] = atom1[i];
        bond_atom2[nbonds] = atom2[i];
        nbonds++;
    }
}

void Topology::add_angls(int n, int *types, int *atom1, int *atom2, int *atom3)
{
    for(int i=0; i<n; i++)
    {
        if(nangls >= maxangls)
        {
            maxangls += 1000;
            angl_types = (int*)realloc(angl_types, sizeof(int) * maxangls);
            angl_atom1 = (int*)realloc(angl_atom1, sizeof(int) * maxangls);
            angl_atom2 = (int*)realloc(angl_atom2, sizeof(int) * maxangls);
            angl_atom3 = (int*)realloc(angl_atom3, sizeof(int) * maxangls);
        }
        
        angl_types[nangls] = types[i];
        angl_atom1[nangls] = atom1[i];
        angl_atom2[nangls] = atom2[i];
        angl_atom3[nangls] = atom3[i];
        nangls++;
    }
}

void Topology::add_dihes(int n, int *types, int *atom1, int *atom2, int *atom3, int *atom4)
{
    for(int i=0; i<n; i++)
    {
        if(ndihes >= maxdihes)
        {
            maxdihes += 1000;
            dihe_types = (int*)realloc(dihe_types, sizeof(int) * maxdihes);
            dihe_atom1 = (int*)realloc(dihe_atom1, sizeof(int) * maxdihes);
            dihe_atom2 = (int*)realloc(dihe_atom2, sizeof(int) * maxdihes);
            dihe_atom3 = (int*)realloc(dihe_atom3, sizeof(int) * maxdihes);
            dihe_atom4 = (int*)realloc(dihe_atom4, sizeof(int) * maxdihes);
        }
        
        dihe_types[ndihes] = types[i];
        dihe_atom1[ndihes] = atom1[i];
        dihe_atom2[ndihes] = atom2[i];
        dihe_atom3[ndihes] = atom3[i];
        dihe_atom4[ndihes] = atom4[i];
        ndihes++;
    }
}

#define MAXSPECIAL 20

void Topology::build_special(bool bonds, bool angles, bool diheds)
{
    nspecials = new int [natoms];
    special_pairs = new int* [natoms];
    int *buf = new int [natoms * MAXSPECIAL];
    
    for(int i=0; i<natoms; i++) 
    {
        nspecials[i] = 0;
        special_pairs[i] = buf + MAXSPECIAL * i;
    }
    
    if(bonds)
    {
        for(int i=0; i<nbonds; i++)
        {
            int ia = bond_atom1[i], ib = bond_atom2[i];
            special_pairs[ia][nspecials[ia]++] = ib;
            special_pairs[ib][nspecials[ib]++] = ia;
        }
    }
    
    if(angles)
    {
        for(int i=0; i<nangls; i++)
        {
            int ia = angl_atom1[i], ib = angl_atom3[i];
            special_pairs[ia][nspecials[ia]++] = ib;
            special_pairs[ib][nspecials[ib]++] = ia;
        }
    }
    
    if(diheds)
    {
        for(int i=0; i<ndihes; i++)
        {
            int ia = dihe_atom1[i], ib = dihe_atom4[i];
            special_pairs[ia][nspecials[ia]++] = ib;
            special_pairs[ib][nspecials[ib]++] = ia;
        }
    }
}





