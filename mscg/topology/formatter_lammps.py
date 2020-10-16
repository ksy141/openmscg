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

def read(rows, **kwargs):
    
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

def write(top, **kwargs):
    
    box = kwargs['box']
    masses = kwargs['masses']
    x = kwargs['x']
    
    with open(top.system_name + '.data', 'w') as f:
        f.write("%s\n\n" % (top.system_name))
        
        for t in Topology.all_types:
            n = len(top._names[t])
            if n>0: f.write("%d %s types\n" % (n, t))
        
        f.write("\n")
        
        for t in Topology.all_types:
            n = len(top._types[t])
            if n>0: f.write("%d %ss\n" % (n, t))
            
        
        f.write("\n%f %f xlo xhi" % (box[0][0], box[0][1]))
        f.write("\n%f %f ylo yhi" % (box[1][0], box[1][1]))
        f.write("\n%f %f zlo zhi" % (box[2][0], box[2][1]))
        f.write("\n")
                
        f.write("\nMasses\n\n")
        
        for i in range(top.ntype_atom):
            f.write("%d %f # %s\n" % (i+1, masses[i], top.names_atom[i]))

        f.write("\nAtoms\n\n")
        types_atom = top.types_atom
        
        for i in range(top.n_atom):
            f.write("%d %s %d %s %f %f %f\n" % (i+1, 
                str(kwargs['molecule'][i]) if 'molecule' in kwargs else '', 
                types_atom[i] + 1,
                str(kwargs['charge'][i]) if 'charge' in kwargs else '', 
                x[i][0], x[i][1], x[i][2]))
        
        for t in Topology.bonded_types:
            n = len(top._types[t])
            
            if n>0:
                f.write("\n%ss\n\n" % (t.capitalize()))
                bonded_atoms = top._bonded_atoms[t]
                
                for i in range(n):
                    atom_index = ' '.join([str(bonded_atoms[d, i] + 1) for d in range(Topology.bonded_natoms[t])])
                    f.write("%d %d %s\n" % (i+1, top._types[t][i]+1, atom_index))

    
    
    