import os
from . import Topology

def read_moltype(rows, ir, systop):
    top = systop.copy()
    w = rows[ir]
    ir += 1
    
    # Read moltype header
    
    if w[0]!='mol':
        raise Exception('Expecting key [mol] at line %d' % (ir))
    
    natoms = int(w[1])
    bondtype = int(w[2])
    
    # Read sitetypes
    
    if rows[ir][0] != "sitetypes":
        raise Exception('Expecting key [sitetypes] at line %d: %s' % (ir, rows[ir]))
    
    ir += 1
        
    top.add_atoms([row[0] for row in rows[ir:ir+natoms]])
    ir += natoms
    
    # Read bonds
    
    w = rows[ir]
    if w[0] != "bonds":
        raise Exception('Expecting key [bonds] at line %d: %s' % (ir, rows[ir]))
    
    ir += 1
    
    nb = int(w[1])
    bonds = [[], []]
    
    for row in rows[ir:ir+nb]:
        bonds[0].append(int(row[0])-1)
        bonds[1].append(int(row[1])-1)
    
    top.add_bondings('bond', ['?'] * nb, bonds, True)
    ir += nb
    
    # Parse bonding information
    
    if bondtype>1:
        top.gen_angles(autotype=True)
        
        if bondtype>2:
            top.gen_dihedrals(autotype=True)
    
    elif bondtype == -1:
        # add explicit angle list
        if rows[ir][0] != "angles":
            raise Exception('Expecting key [angles] at line %d: %s' % (ir, rows[ir]))
        
        na = int(rows[ir][1])
        fa = int(rows[ir][2])
        a1, a2, a3 = [], [], []
        ir += 1
        
        for row in rows[ir:ir+na]:
            a1.append(int(row[0])-1)
            a2.append(int(row[1])-1)
            a3.append(int(row[2])-1)
        
        if fa == 1:
            top.add_bondings('angle', ['?'] * na, [a1, a2, a3], True)
        else:
            top.add_bondings('angle', ['?'] * na, [a2, a1, a3], True)
        
        ir += na
        
        # add explicit dihedral list
        if rows[ir][0] != "dihedrals":
            raise Exception('Expecting key [dihedrals] at line %d: %s' % (ir, rows[ir]))
        
        nd = int(rows[ir][1])
        fd = int(rows[ir][2])
        d1, d2, d3, d4 = [], [], [], []
        ir += 1
        
        for row in rows[ir:ir+nd]:
            d1.append(int(row[0])-1)
            d2.append(int(row[1])-1)
            d3.append(int(row[2])-1)
            d4.append(int(row[3])-1)
        
        if fd == 1:
            top.add_bondings('dihedral', ['?'] * nd, [d1, d2, d3, d4], True)
        else:
            top.add_bondings('dihedral', ['?'] * nd, [a3, a1, a2, a4], True)
        
        ir += nd    
    
    return ir, top

def inspect(rows):
    keys = [row[0] for row in rows if len(row)>0]
    
    if 'cgsites' in keys and 'cgtypes' in keys:
        return True
    
    return False
    
def read(rows):
    top = Topology()
    ir = 0
    
    # parse file
    
    atoms = []
    bonds = []
    angls = []
    dihes = []
    
    while ir < len(rows):
        w = rows[ir]
        
        if rows[ir] == [] or w[0] in ['cgsites']:
            ir += 1
        
        elif w[0] == 'cgtypes':
            n = int(w[1])
            top.add_names('atom', [row[0] for row in rows[ir+1:ir+n+1]])
            ir += n + 1
        
        elif w[0] == 'moltypes':
            nmt = int(w[1])
            moltypes = [None] * nmt
            ir += 1
            
            for i in range(int(w[1])):
                ir, moltypes[i] = read_moltype(rows, ir, top)
        
        elif w[0] == 'system':
            nm = int(w[1])
            ir += 1
            start = 0
            
            for im in range(nm):
                w = [int(s) for s in rows[ir]]
                ir += 1
                
                mt = moltypes[w[0]-1]
                top += mt * w[1]
            
        else:
            raise Exception('Error parsing line %d: %s' % (ir+1, rows[ir]))
    
    return top














