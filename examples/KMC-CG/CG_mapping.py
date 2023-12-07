import yaml, os, MDAnalysis as mda, numpy as np, pandas as pd
from mscg.cli import cgmap, cghenm
import argparse

# reading topology and crd from mdanalysis
# these variables are changeable according to your needs

parser = argparse.ArgumentParser(description="Generates CG mapped AA trajectory, CG-HENM force field and prepares the LAMMPS input file.")
parser.add_argument("-top", "--top", type=str, help = 'Topology file of the AA trajectory')
parser.add_argument("-crd", "--crd", type=str, help="**ALIGNED** trajectory file of the AA trajectory")
parser.add_argument("-cg", "--cgmap", "--map", type=str, help="CG mapping result file.")
parser.add_argument("-o", "--out", type=str, help="User-set prefix for output files.")
args = parser.parse_args()

"""topo_file = "traj/G-ATP-Num-12/protMGnuc.psf"
crd_file = "traj/G-ATP-Num-12/protMGnuc_aligned.dcd"
cgfile = "traj/G-ATP-Num-12/best_results/result_310.txt"
outfile = "G-ATP-12_310"""

topo_file = args.top
crd_file = args.crd
cgfile = args.cgmap
outfile = args.out

system_info = mda.Universe(topo_file, crd_file)
protein_info = system_info.atoms
residue_info = protein_info.split('residue')

# reading mapping labels from file

with open(cgfile) as f:
    aa = f.readlines()
    mapping_labels = []
    for line in aa:
        try:
            lst = eval(line)
            lst = list(np.array(lst))
            lst[0]
            mapping_labels.append(lst)
        except:
            pass

# Calculating sum of mass for each cg bead

masses = []
for cg_bead in mapping_labels:
    sum_of_mass = 0
    for residue in cg_bead:
        sum_of_mass += np.sum(residue_info[residue].masses)
    masses.append(sum_of_mass)

# loading yaml file

num_CG_site = len(mapping_labels)
offset = len(protein_info)
map_output = {"site-types":{}, "system":[]}
for CG_site_index in range(num_CG_site):
    map_output["site-types"]["CG{}".format(CG_site_index+1)] = {}
    map_output["site-types"]["CG{}".format(CG_site_index+1)]['index'] = []
    map_output["site-types"]["CG{}".format(CG_site_index+1)]['x-weight'] = []
    map_output["site-types"]["CG{}".format(CG_site_index+1)]['f-weight'] = []
    '''if CG_site_index >= num_CG_site - 2: #PA PB PG
        for i in mapping_labels[CG_site_index]:
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['index'].append(int(residue_info[i].select_atoms("name PA or name PB or name PG").indices[0]))
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['x-weight'].append(30.9738)
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['f-weight'].append(1.0)
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['index'].append(int(residue_info[i].select_atoms("name PA or name PB or name PG").indices[1]))
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['x-weight'].append(30.9738)
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['f-weight'].append(1.0)
    else: #CA'''
    ###start
    if 1:
    ###end
        for i in mapping_labels[CG_site_index]:
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['index'].append(int(residue_info[i].select_atoms("name CA").indices[0]))
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['x-weight'].append(12.011)
            map_output["site-types"]["CG{}".format(CG_site_index+1)]['f-weight'].append(1.0)

map_output["system"].append({})
map_output["system"][0]["anchor"] = 0
map_output["system"][0]["repeat"] = 1
map_output["system"][0]["offset"] = offset
map_output["system"][0]["sites"] = []

for CG_site_index in range(num_CG_site):
    map_output["system"][0]["sites"].append(["CG{}".format(CG_site_index+1), 0])

# writing yaml file

with open("%s.yaml" % outfile, 'w') as ff:
    yaml.dump(map_output, ff, default_flow_style=None)

# cgmapping of the all_atom trajectory

cgmap.main(map = "%s.yaml" % outfile, traj = crd_file, out="%s.lammpstrj" % outfile)

# cghenm modelling

cghenm.main(traj = "%s.lammpstrj" % outfile, temp = 310.0, cut=75.0, alpha=0.5, save="%s_cghenm" % outfile)

# Reading initial CG coordinates from lammpstrj file

assert "%s_cghenm.txt" % outfile in os.listdir(".")
with open( "%s.lammpstrj" % outfile) as trj:
    while True:
        a = trj.readline()
        if "BOX BOUNDS" in a: break
    box = []
    i=0
    while i<3: 
        box.append(list(np.asarray(trj.readline().strip().split(), dtype=float)))
        i+=1
    while True:
        a = trj.readline()
        if "ATOMS" in a: 
            break
    atom_ids, atom_types, initial_coords = [], [], []
    while True:
        a = trj.readline()
        if "ITEM" in a: break
        info = a.strip().split()
        atom_ids.append(int(info[0]))
        atom_types.append(int(info[1]))
        initial_coords.append(info[-3:])

# Generating lammps data file

with open("%s_cghenm.txt" % outfile) as henm:
    aa = henm.readlines()
    info = [line.strip().split() for line in aa[1:]]
    info = np.asarray(info, dtype=float)
    header_list = aa[0].strip().split("  ")
    header_list = [item.strip() for item in header_list]
    atom_pairs = np.asarray(info[:,:2], dtype=int)
    coeffs = np.asarray(info[:,2:4], dtype=float)

with open("%s.lammpsdat" % outfile, 'w') as f:
    print("File generated by CG cluster program\n", file = f)
    print(len(masses), "atoms", file=f)
    print("%d bonds\n0 angles\n0 dihedrals\n" % len(atom_pairs), file = f)
    print(len(masses), "atom types", file=f)
    print("%d bond types\n0 angle types\n0 dihedral types\n" % len(atom_pairs), file = f)
    print("%.6f %.6f" % (box[0][0], box[0][1]), "xlo xhi", file=f)
    print("%.6f %.6f" % (box[1][0], box[1][1]), "ylo yhi", file=f)
    print("%.6f %.6f" % (box[2][0], box[2][1]), "zlo zhi", file=f)

    print("\nMasses\n", file=f)
    for i in range(len(masses)):
        print(i+1, "%.3f" % masses[i], file=f)

    print("\nAtoms # full\n", file=f)
    for i in range(len(masses)):
        print("%d 1 %d 0 %s" % (atom_ids[i], atom_types[i], ' '.join(initial_coords[i])), file=f)
        # 0 stands for the charge

    print("\nBonds\n", file=f)
    for i in range(len(atom_pairs)):
        print(i+1, i+1, "%d %d" % (atom_pairs[i][0], atom_pairs[i][1]), file=f)

    print("\nBond Coeffs\n", file=f)
    for i in range(len(atom_pairs)):
        print(i+1, "%.6f %.6f" % (coeffs[i][1], coeffs[i][0]), file=f) # lammps default arrangement: R0 K


with open("%s.in" % outfile, 'w') as infile:
    print("""########### CG simulation ###########
clear

variable TotalTime equal 100000000

units \t real
dimension \t 3
atom_style \t full
bond_style \t harmonic
boundary \t p p p

read_data \t {0}.lammpsdat

timestep \t 10
run_style \t verlet
run \t 0

reset_timestep \t 0

velocity \t all create 310.0 2022416 mom yes rot yes
# "mom yes" means initial momentum is set to zero.
# "rot yes" means initial angular momentum is set to zero.
fix \t 1 all langevin 310.0 310.0 1000 2022327
fix \t 2 all nve

thermo_style \t custom step temp press ke pe evdwl etotal
thermo \t 20000
# outputs the thermo information every 20000 steps.

dump this_is_just_a_name all custom 20000 {0}_outtraj.lammpstrj id type x y z

run \t ${{TotalTime}}
write_data {0}.dat
""".format(outfile), file=infile)
