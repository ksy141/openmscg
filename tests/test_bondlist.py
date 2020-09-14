import pytest
from mscg import *

def test_glu(datafile):
    top = Topology.read_file('tests/data/lammps_glu.data')
    trj = Trajectory('tests/data/lammps_glu.lammpstrj', fmt='lammpstrj')
    trj.read_frame()
    
    blist = BondList(
        top.types_bond, top.bond_atoms, 
        top.types_angle, top.angle_atoms, 
        top.types_dihedral, top.dihedral_atoms)
    
    blist.build(trj.box, trj.x)
    
    assert abs(blist.get_scalar('bond').mean() - 1.1142256) < 0.001
    assert abs(blist.get_scalar('angle').mean() - 1.9674782) < 0.001
