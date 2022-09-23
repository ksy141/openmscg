'''Solvation free enery mapping

Description
-----------

The ``cgsolv`` command is used to ...

https://amberhub.chpc.utah.edu/gist/

 0 • index – A unique, sequential integer assigned to each voxel
 1 • xcoord – x coordinate of the center of the voxel (Å)
 2 • ycoord – y coordinate of the center of the voxel (Å)
 3 • zcoord – z coordinate of the center of the voxel (Å)
 4 • population – Number of water molecule, ni, found in the voxel over the entire simulation. A water molecule is deemed to populate a voxel if its oxygen coordinates are inside the voxel. The expectation value of this quantity increases in proportion to the length of the simulation.
 5 • g_O – Number density of oxygen centers found in the voxel, in units of the bulk density (rdval). Thus, the expectation value of g_O for a neat water system is unity.
 6 • g_H – Number density of hydrogen centers found in the voxel in units of the reference bulk density (2rdval). Thus, the expectation value of g_H for a neat water system would be unity.
 7 • dTStrans-dens – First order translational entropy density (kcal/mole/Å3), referenced to the translational entropy of bulk water, based on the value rdval.
 8 • dTStrans-norm – First order translational entropy per water molecule (kcal/mole/molecule), referenced to the translational entropy of bulk water, based on the value rdval. The quantity dTStrans-norm equals dTStrans-dens divided by the number density of the voxel.
 9 • dTSorient-dens – First order orientational entropy density (kcal/mole/Å3), referenced to bulk solvent (see below).
10 • dTSorient-norm – First order orientational entropy per water molecule (kcal/mole/water), referenced to bulk solvent (see below). This quantity equals dTSorient-dens divided by the number density of the voxel.
11 • Esw-dens – Mean solute-water interaction energy density (kcal/mole/Å3). This is the interaction of the solvent in a given voxel with the entire solute. Both Lennard-Jones and electrostatic interactions are computed without any cutoff, within the minimum image convention but without Ewald summation. This quantity is referenced to bulk, in the trivial sense that the solute-solvent interaction energy is zero in bulk.
12 • Esw-norm – Mean solute-water interaction energy per water molecule. This equals Esw-dens divided by the number density of the voxel (kcal/mole/molecule).
13 • Eww-dens – Mean water-water interaction energy density, scaled by ½ to prevent double-counting, and not referenced to the corresponding bulk value of this quantity (see below). This quantity is one half of the mean interaction energy of the water in a given voxel with all other waters in the system, both on and off the GIST grid, divided by the volume of the voxel (kcal/mole/Å3). Again, both Lennard-Jones and electrostatic interactions are computed without any cutoff, within the minimum image convention.
14 • Eww-norm – Mean water-water interaction energy, normalized to the mean number of water molecules in the voxel (kcal/mole/water). See prior column definition for details.
15 • Dipole_x-dens – x-component of the mean water dipole moment density (Debye/Å3).
16 • Dipole_y-dens – y-component of the mean water dipole moment density (Debye/Å3).
17 • Dipole_z-dens – z-component of the mean water dipole moment density (Debye/Å3).
18 • Dipole-dens – Magnitude of mean dipole moment (polarization) (Debye/Å3).
19 • Neighbor-dens – Mean number of waters neighboring the water molecules

Usage
-----

Syntax of running ``cgib`` command ::

    usage: cgib [-h] [-v] --top file [--names] [--traj file[,args]] [--cut]
            [--temp] [--pair types,args] [--bond types,args]
            [--angle types,args] [--plot U or N]

    General arguments:
      -h, --help          show this help message and exit
      -v , --verbose      screen verbose level (default: 0)

'''

from mscg import *
from mscg.pdb import PDB
from mscg.mass import Masses

import pandas as pd

def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'CGSOLV please read ' + doc_root + 'commands/cgsolv.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--pdb",  metavar='file', type=str, help="protein PDB file", required=True)
    group.add_argument("--map",  metavar='file', type=str, help="Atom-to-CG mapping file", required=True)
    group.add_argument("--gist",  metavar='file', type=str, help="GIST output file", required=True)
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--pmin", metavar='', type=float, default=1.5, help="threshold for population numbers")
    group.add_argument("--eww", metavar='', type=float, default=-0.3704, help="bulk water solvation energy")
    group.add_argument("--rmax", metavar='', type=float, default=30.0, help="cut-off for solvation energy calculations")
    group.add_argument("--bin", metavar='', type=float, default=0.1, help="width of bins")
    
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
        
    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)
    
    # load pdb and build CG coordinates
    screen.info("Read solute structure data ...")
    pdb = PDB.read(args.pdb)
    
    with open(args.map, 'r') as f:
        seq = f.read().strip().split("\n")
    
    pdb.atoms['CGID'] = pd.Series(seq)
    w = Masses.get(list(pdb.atoms['ELEMENT']))
    mapper = Mapper.build_from_sequence(seq, w, [1.0] * len(w))
        
    x_solute = pdb.get_coords()
    x_cg = mapper.process(np.array([0.0]*3), x_solute, None)[0]
    cg = pd.DataFrame(x_cg, columns=['XCG','YCG','ZCG'])
    cg['CGID'] = pd.Series([_[0] for _ in mapper.sites])    
    
    # load gist output
    screen.info("Read and process GIST solvation grid ...")
    gdata = np.loadtxt(args.gist).T
    gdata[0, :] = np.arange(0, gdata.shape[1])
    gdata = gdata[:, gdata[4] > args.pmin]
    esolv = gdata[13]*0.5 + gdata[15] - gdata[7] - gdata[9] - args.eww
    gdata = gdata[0:5]
    gdata[4, :] = esolv
    
    # solvent grid assignment
    def search_grid(gdata, x_solute, i):
        if i%100==0: print(i)
        
        dr_v = gdata[1:4].T - x_solute[i]
        dr_s = np.sqrt(np.sum(np.square(dr_v), axis=1))
        
        sel = dr_s<args.rmax
        gdata = gdata[:3, sel]
        gdata[1, :] = dr_s[sel]
        gdata[2, :] = np.full(gdata.shape[1], i)
        return gdata.T
    
    # build dG for each CG siite    
    assoc = np.vstack([search_grid(gdata, x_solute, i) for i in range(x_solute.shape[0])])
    assoc = pd.DataFrame(assoc, columns=['ID','R','SID'])
    assoc = assoc.sort_values('R').groupby('ID').first().reset_index()
    
    assoc = assoc.merge(pd.DataFrame(gdata[:5].T, columns=['ID','X','Y','Z','W']), on='ID', how='left') \
        .merge(pd.DataFrame({'SID':np.arange(x_solute.shape[0]),'CGID':seq}), on='SID', how='left') \
        .merge(cg, on='CGID', how='left')
    
    assoc = assoc[assoc.R>1.5]
    assoc['RCG'] = np.sqrt((assoc['X'] - assoc['XCG'])**2 + (assoc['Y'] - assoc['YCG'])**2 + (assoc['Z'] - assoc['ZCG'])**2)
    assoc['TABLE'] = np.round(assoc.RCG / args.bin) * args.bin
    assoc.CGID = assoc.CGID.astype(int)
    print(assoc)
    
    tables = assoc.groupby(['CGID','TABLE']).agg({'W':'mean', 'TABLE':'count'})
    tables['V1'] = (tables['TABLE']/8) / (np.pi*((tables['W']+0.5).pow(3.0)-(tables['W']-0.5).pow(3.0))/9)
    tables['V2'] = 1.0
    tables['W'] *= tables[['V1', 'V2']].min(axis=1)
    print(tables)
    
    idx = pd.MultiIndex.from_product([
        assoc.CGID.unique(), list(np.arange(0.0, args.rmax + args.bin, args.bin))], 
        names=["CGID", "TABLE"])
    
    dG = pd.DataFrame(index = idx)
    dG['W'] = tables['W']
    dG = dG.W.fillna(0.0).unstack().fillna(0.0)
    print(dG)
    np.savetxt("dG.dat", dG.values)
    
    
if __name__ == '__main__':
    main()
