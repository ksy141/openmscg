import os

from . import Topology

def inspect(rows):
    has_masses, has_atom_types = False, False
    
    for row in rows:
        if 'Masses' in row:
            has_masses = True
        
        if "types" in row:
            has_atom_types = True
        
        if has_masses and has_atom_types:
            return True
    
    return False

def load(rows):
    
    def find_section(rows, header):
        try:
            section_start = rows.index(header) + 2
        except:
            return []
        
        section_end = section_start
        
        while section_end<len(rows) and rows[section_end]!=[]:
            section_end += 1
        
        return range(section_start, section_end)
    
    top = Topology()
    top.add_names('atom', [rows[i][0] for i in find_section(rows, ['Masses'])])
    top.add_atoms([rows[i][2] for i in find_section(rows, ['Atoms'])])
    
    for t in Topology.bonded_types:
        atoms = [[] for _ in range(Topology.bonded_natoms[t])]
        
        for i in find_section(rows, [t.capitalize() + 's']):
            for j in range(Topology.bonded_natoms[t]):
                atoms[j].append(int(rows[i][2+j])-1)
        
        top.add_bondings(t, ['?'] * len(atoms[0]), atoms, autotype=True)
    
    return top
