import os
from .topology import *

def read_moltype(rows, ir):
    
    w = rows[ir].split(" ")
    ir += 1
    
    # Read moltype header
    
    if w[0]!='mol':
        raise Exception('Expecting key [mol] at line %d' % (ir))
    
    natoms = int(w[1])
    bondtype = int(w[2])
    
    # Read sitetypes
    
    if rows[ir].strip() != "sitetypes":
        raise Exception('Expecting key [sitetypes] at line %d: %s' % (ir, rows[ir]))
    
    ir += 1
    
    m = {}
    m['atoms'] = [int(row)-1 for row in rows[ir:ir+natoms]]
    ir += natoms
    
    # Read bonds
    
    w = rows[ir].split(" ")
    if w[0] != "bonds":
        raise Exception('Expecting key [bonds] at line %d: %s' % (ir, rows[ir]))
    
    ir += 1
    
    nb = int(w[1])
    m['bonds'] = [[int(s)-1 for s in row.split(" ")] for row in rows[ir:ir+nb]]
    ir += nb
    
    # Parse bonding information
    
    if bondtype>1:
        m['angls'] = generate_angles(len(m['atoms']), m['bonds'])
        
        if bondtype>2:
            m['dihes'] = generate_dihedrals(len(m['atoms']), m['bonds'])
    
    elif bondtype == -1:
        # Manual angles and dihedrals here
        pass
    
    return ir, m
    
    
    

def build_top_from_cgtop(filename):

    # read file content

    if not os.path.isfile(filename):
        raise Exception('File is not found: ' + filename)

    with open(filename, 'r') as f:
        rows = f.read().split("\n")
    
    # parse file
    
    atoms = []
    bonds = []
    angls = []
    dihes = []
            
    top = Topology()
    ir = 0
    
    while ir < len(rows):
        w = [s.strip() for s in rows[ir].split(" ")]
        
        if rows[ir] == "" or w[0] in ['cgsites']:
            ir += 1
        
        elif w[0] == 'cgtypes':
            n = int(w[1])
            top.set_names(rows[ir+1:ir+n+1])
            ir += n + 1
        
        elif w[0] == 'moltypes':
            nmt = int(w[1])
            moltypes = [None] * nmt
            ir += 1
            
            for i in range(int(w[1])):
                ir, moltypes[i] = read_moltype(rows, ir)
        
        elif w[0] == 'system':
            nm = int(w[1])
            ir += 1
            start = 0
            
            for im in range(nm):
                w = [int(s) for s in rows[ir].split(" ")]
                ir += 1
                
                mt = moltypes[w[0]-1]
                mols = w[1]
                
                for _ in range(mols):
                    atoms.extend([i for i in mt['atoms']])
                    bonds.extend([[b[0]+start, b[1]+start] for b in mt['bonds']])
                    angls.extend([[a[0]+start, a[1]+start, a[2]+start] for a in mt['angls']])
                    dihes.extend([[d[0]+start, d[1]+start, d[2]+start, d[3]+start] for d in mt['dihes']])
                    start += len(moltypes[w[0]-1]['atoms'])
            
        else:
            raise Exception('Error parsing line %d: %s' % (ir+1, rows[ir]))
    '''
    print("Atoms:", atoms)
    print("Bonds:", bonds)
    print("Angles:", angls)
    print("Dihedrals:", dihes)
    '''
    # create c++ object
    
    top.create()
    top.add_atoms(atoms)
    top.add_bonds([i[0] for i in bonds], [i[1] for i in bonds])
    top.add_angles([i[0] for i in angls], [i[1] for i in angls], [i[2] for i in angls])
    top.add_dihedrals([i[0] for i in dihes], [i[1] for i in dihes], [i[2] for i in dihes], [i[3] for i in dihes])
    
    return top














