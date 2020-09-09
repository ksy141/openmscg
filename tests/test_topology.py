import pytest
from mscg import Topology

def _build_waterbox(nwater):
    top = Topology()
    
    top.add_names('atom', ['OW', 'HW'])
    top.add_names('bond', 'OW-HW')
    top.add_names('angle', 'HW-OW-HW')
    
    for i in range(nwater):
        n = top.add_atoms(['OW', 'HW', 'HW'])
        top.add_bondings('bond', ['?'] * 2, [[n-3] * 2, [n-2, n-1]], True)
        top.add_bondings('angle', ['?'], [[n-2], [n-3], [n-1]], True)
    
    return top

def _build_methanol(nmol):
    top = Topology()
    
    top.add_names('atom', ['CT', 'HT', 'OH', 'HO'])
        
    for i in range(nmol):
        n = top.add_atoms(['CT'] + ['HT']*3 + ['OH', 'HO'])
        top.add_bondings('bond', ['?'] * 5, [
            [n-6] * 4 + [n-2], 
            [n-5, n-4, n-3, n-2, n-1]
        ], autotype = True)
    
    top.gen_angles(autotype = True)
    top.gen_dihedrals(autotype = True)
    
    return top
    
def test_waterbox():
    nwater = 8
    top = _build_waterbox(nwater)
    
    top.add_names('atom', ['OW'])
    assert len(top.names_atom) == 2
    
    assert top.n_atom == nwater * 3
    assert top.n_bond == nwater * 2
    assert top.n_angle == nwater
    
    links = top.linking_map()
    assert links.shape == (nwater * 3, 3)
    assert (links!=-1).sum() == nwater * 6
    
    assert top.pair_type([0, 1]) == (1, 'OW-HW')
    
def test_autogen():
    nmol = 5
    top = _build_methanol(nmol)
    
    assert top.names_bond == ['CT-HT', 'CT-OH', 'OH-HO']
    assert top.n_bond == nmol * 5
    assert sum(top.types_bond) == nmol * 3
    
    assert top.names_angle == ['HT-CT-HT', 'HT-CT-OH', 'CT-OH-HO']
    assert top.n_angle == nmol * 7
    assert sum(top.types_angle) == nmol * 5
    
    assert top.names_dihedral == ['HT-CT-OH-HO']
    assert top.n_dihedral == nmol * 3
    assert sum(top.types_dihedral) ==0
    
    top.reset_names(['C1', 'H1', 'O', 'H'])
    assert top.names_bond == ['C1-H1', 'C1-O', 'O-H']
    assert top.names_angle == ['H1-C1-H1', 'H1-C1-O', 'C1-O-H']
    assert top.names_dihedral == ['H1-C1-O-H']

def test_overload():
    # copy
    water = _build_waterbox(1)
    top_copy = water.copy()
    assert water.types_atom.sum() == top_copy.types_atom.sum()
    assert water.types_bond.sum() == top_copy.types_bond.sum()
    assert water.bond_atoms.sum() == top_copy.bond_atoms.sum()
    
    # multiply
    top_mul = water * 3
    assert top_mul.types_atom.sum() == water.types_atom.sum() * 3
    assert top_mul.types_bond.sum() == water.types_bond.sum() * 3
    assert top_mul.bond_atoms.sum() == 45
    
    # self-add
    methanol = _build_methanol(1)
    top_copy += methanol
    assert top_copy.names_atom == ['OW', 'HW', 'CT', 'HT', 'OH', 'HO']
    assert top_copy.names_bond == ['OW-HW', 'CT-HT', 'CT-OH', 'OH-HO']
    assert top_copy.types_atom.sum() == 22
    assert top_copy.types_bond.sum() == 8
    assert top_copy.bond_atoms.sum() == 52
    
    # add
    top = methanol * 2 + water
    assert top.names_atom == ['CT', 'HT', 'OH', 'HO', 'OW', 'HW']
    assert top.names_bond == ['CT-HT', 'CT-OH', 'OH-HO', 'OW-HW']
    assert top.types_atom.sum() == 30
    assert top.types_bond.sum() == 12
    assert top.bond_atoms.sum() == 149
    
def test_formatter_lammps(datafile):
    top = Topology.read_file(datafile('lammps_glu.data'))
    assert len(top.names_dihedral) == 39
    assert len(top.types_dihedral) == 60
    assert top.types_angle.sum() == 581
    
def test_formatter_cgtop(datafile):
    top = Topology.read_file(datafile('4_site.top'))
    assert [top.n_atom, top.n_bond, top.n_angle, top.n_dihedral] == [4, 3, 2, 1]
    assert top.names_angle == ['1-1-1']
