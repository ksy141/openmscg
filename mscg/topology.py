#

import os, sys
from .core import cxx_topol as lib

class Topology:

    def __init__(self):
        self.h = None

        self.ntypes_atom, self.ntypes_bond, self.ntypes_angl, self.ntypes_dihe = 0, 0, 0, 0
        self.types_atom, self.types_bond, self.types_angl, self.types_dihe = [], [], [], []

        self.natoms, self.nbonds, self.nangls, self.ndihes = 0, 0, 0, 0
        self.atom_types = []
        self.bond_types, self.bond_atom1, self.bond_atom2 = [], [], []
        self.angl_types, self.angl_atom1, self.angl_atom2, self.angl_atom3 = [], [], [], []
        self.dihe_types, self.dihe_atom1, self.dihe_atom2, self.dihe_atom3, self.dihe_atom4 = [], [], [], [], []



    def __del__(self):
        if self.h is not None:
            lib.destroy(self.h)

    def create(self):
        if self.h is not None:
            lib.destroy(self.h)

        self.h = lib.create(self.ntypes_atom)

    def set_names(self, atom_names):
        self.types_atom = atom_names
        self.ntypes_atom = len(atom_names)

    def add_atoms(self, atom_types):
        if type(atom_types[0]) == str:
            atom_types = [self.get_atom_type(i) for i in atom_types]

        self.atom_types.extend(atom_types)
        lib.add_atoms(self.h, atom_types)
        self.natoms = len(self.atom_types)

    def add_bonds(self, bond_atom1, bond_atom2):
        types = self.types_atom
        bond_types = []

        for i in range(len(bond_atom1)):
            a1 = self.atom_types[bond_atom1[i]]
            a2 = self.atom_types[bond_atom2[i]]
            bond = self.get_bond_name(types[a1], types[a2])

            if bond not in self.types_bond:
                self.types_bond.append(bond)
                bond_types.append(len(self.types_bond)-1)
            else:
                bond_types.append(self.get_bond_type(bond))

        self.bond_types.extend(bond_types)
        self.bond_atom1.extend(bond_atom1)
        self.bond_atom2.extend(bond_atom2)
        self.ntypes_bond = len(self.types_bond)
        lib.add_bonds(self.h, bond_types, bond_atom1, bond_atom2)

    def add_angles(self, atom1, atom2, atom3):
        types = self.types_atom
        angl_types = []

        for i in range(len(atom1)):
            a1 = self.atom_types[atom1[i]]
            a2 = self.atom_types[atom2[i]]
            a3 = self.atom_types[atom3[i]]
            angle = self.get_angle_name(types[a1], types[a2], types[a3])

            if angle not in self.types_angl:
                self.types_angl.append(angle)
                angl_types.append(len(self.types_angl)-1)
            else:
                angl_types.append(self.get_angle_type(angle))

        self.angl_types.extend(angl_types)
        self.angl_atom1.extend(atom1)
        self.angl_atom2.extend(atom2)
        self.angl_atom3.extend(atom3)
        self.ntypes_angl = len(self.types_angl)
        lib.add_angles(self.h, angl_types, atom1, atom2, atom3)

    def add_dihedrals(self, atom1, atom2, atom3, atom4):
        types = self.types_atom
        dihe_types = []

        for i in range(len(atom1)):
            a1 = self.atom_types[atom1[i]]
            a2 = self.atom_types[atom2[i]]
            a3 = self.atom_types[atom3[i]]
            a4 = self.atom_types[atom4[i]]
            di = self.get_dihedral_name(types[a1], types[a2], types[a3], types[a4])

            if di not in self.types_dihe:
                self.types_dihe.append(di)
                dihe_types.append(len(self.types_dihe)-1)
            else:
                dihe_types.append(self.get_dihedral_type(di))

        self.dihe_types.extend(dihe_types)
        self.dihe_atom1.extend(atom1)
        self.dihe_atom2.extend(atom2)
        self.dihe_atom3.extend(atom3)
        self.dihe_atom4.extend(atom4)
        self.ntypes_dihe = len(self.types_dihe)
        lib.add_dihedrals(self.h, dihe_types, atom1, atom2, atom3, atom4)

    def get_atom_type(self, name):
        return self.types_atom.index(name)

    def get_bond_type(self, n1, n2=None):
        if n2 is None:
            name = n1
        else:
            name = self.get_bond_name(n1, n2)
        return self.types_bond.index(name)

    def get_bond_name(self, n1, n2):
        return  n1 + '-' + n2 if n1 < n2 else n2 + '-' + n1

    def get_angle_type(self, n1, n2=None, n3=None):
        if n2 is None:
            name = n1
        else:
            name = self.get_angle_name(n1, n2, n3)
        return self.types_angl.index(name)

    def get_angle_name(self, n1, n2, n3):
        seq = [n1, n2, n3] if n1<n3 else [n3, n2, n1]
        return '-'.join(seq)

    def get_dihedral_type(self, n1, n2=None, n3=None, n4=None):
        if n2 is None:
            name = n1
        else:
            name = self.get_dihedral_name(n1, n2, n3, n4)
        return self.types_dihe.index(name)

    def get_dihedral_name(self, n1, n2, n3, n4):
        seq = [n1, n2, n3, n4] if n2<n3 else [n4, n3, n2, n1]
        return '-'.join(seq)

    def get_pair_type(self, nameA, nameB):
        a = self.get_atom_type(nameA)
        b = self.get_atom_type(nameB)
        return a * self.ntypes_atom + b if a<b else b * self.ntypes_atom + a

    def build_special(self, use_bond = False, use_angle = False, use_dihed = False):
        if self.h is None:
            raise Exception('Topology does not exist.')

        lib.build_special(self.h, use_bond, use_angle, use_dihed)
    
    def reset_names(self, names):
        
        if len(names) != self.ntypes_atom:
            raise ValueException("inconsistent number of types in reset_names()")
        
        m = {}
        
        for i in range(self.ntypes_atom):
            m[self.types_atom[i]] = names[i]
        
        self.types_atom = names
        
        t = [one.split('-') for one in self.types_bond]
        self.types_bond = ['-'.join([m[tt[0]], m[tt[1]]]) for tt in t]
        
        t = [one.split('-') for one in self.types_angl]
        self.types_angl = ['-'.join([m[tt[0]], m[tt[1]], m[tt[2]]]) for tt in t]
        
        t = [one.split('-') for one in self.types_dihe]
        self.types_dihe = ['-'.join([m[tt[0]], m[tt[1]], m[tt[2]], m[tt[3]]]) for tt in t]



def generate_angles(natoms, bonds):
    links = [[] for _ in range(natoms)]

    for b in bonds:
        links[b[0]].append(b[1])
        links[b[1]].append(b[0])

    angls = []

    for i in range(natoms):
        link = links[i]
        n = len(link)

        if n<2: continue

        for j in range(n-1):
            for k in range(1, n):
                angls.append([link[j], i, link[k]])

    return angls



def generate_dihedrals(natoms, bonds):
    links = [[] for _ in range(natoms)]

    for b in bonds:
        links[b[0]].append(b[1])
        links[b[1]].append(b[0])

    dihes = []

    for bond in bonds:
        i, ni = bond[0], len(links[bond[0]])
        j, nj = bond[1], len(links[bond[1]])
        if ni<2 or nj<2: continue

        for k in range(ni):
            for l in range(nj):
                if links[i][k] != j and links[j][l] != i:
                    dihes.append([links[i][k], i, j, links[j][l]])

    return dihes

from .top_lammps import *
from .top_cg import*

def build_topology(file):
    
    with open(file, "r") as f:
        text = f.read()
    
    segs = file.split(".")
    suffix = segs[-1] if len(segs)>1 else ""
    
    if ("atom types" in text) and ("Masses" in text):
        return build_top_from_lammps(file)
    elif ("cgsites " in text) and ("cgtypes " in text):
        return build_top_from_cgtop(file)
    else:
        raise Exception("Unknown topology type: " + fmt)
