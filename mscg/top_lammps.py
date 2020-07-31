
from .topology import *

def build_top_from_lammps(filename):

    # read file content

    if not os.path.isfile(filename):
        raise Exception('File is not found: ' + filename)

    with open(filename, 'r') as f:
        rows = [row.split('#')[0].strip().split() for row in f.read().split("\n")]
    
    ir = 2
    
    # read masses (atom types)
    
    masses = []
    
    while ir < len(rows):
        if(len(rows[ir]) > 0):
            if(rows[ir][0] == 'Masses'):
                break
        ir += 1

    ir += 2

    while ir < len(rows) and rows[ir]!=[]:
        row = rows[ir]
        if len(row)<2: raise Exception("Wrong format for [Masses] block at line %d." % (ir))
        masses.append(row[1])
        ir += 1
    
    # read atoms

    atoms = []

    while ir < len(rows):
        if(len(rows[ir]) > 0):
            if(rows[ir][0] == 'Atoms'):
                break
        ir += 1

    ir += 2

    while ir < len(rows) and rows[ir]!=[]:
        row = rows[ir]
        if len(row)<3: raise Exception("Wrong format for [Atoms] block at line %d." % (ir))
        atoms.append(int(row[2])-1)
        ir += 1
    
    # read bonds

    bond_atom1, bond_atom2 = [], []

    while ir < len(rows) and rows[ir]!= ['Bonds']:
        ir += 1

    ir += 2

    while ir < len(rows) and rows[ir]!=[]:
        row = rows[ir]
        if len(row)!=4: raise Exception("Wrong format for [Bonds] block at line %d." % (ir))
        bond_atom1.append(int(row[2])-1)
        bond_atom2.append(int(row[3])-1)
        ir += 1
    
    # read angles

    #

    # read dihedrals

    #

    # create c++ handle
    
    top = Topology()
    top.set_names(masses)
    top.create()
    top.add_atoms(atoms)
    top.add_bonds(bond_atom1, bond_atom2)
    
    if top.ntypes_atom != len(top.types_atom):
        raise Exception("Count of atom names is not consistent.")

    if top.atom_types == []: 
        raise Exception("No atom is defined.")
    
    return top
