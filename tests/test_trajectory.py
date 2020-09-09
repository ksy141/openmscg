import pytest
import numpy as np
from mscg import Trajectory, TrajReader
    
def test_trr(datafile):
    trj = Trajectory(datafile('methanol_1728_cg.trr'), fmt='trr')
    assert trj.natoms == 1728
    assert [trj.has_attr(x) for x in 'tqxvf'] == [False, False, True, False, True]
    
    trj.read_frame()
    assert trj.x.shape == (1728, 3)
    assert abs(trj.x.sum() - 127829.11) < 0.01

def test_lmp(datafile):
    trj = Trajectory(datafile('unary_lj_fluid.lammpstrj'), fmt='lammpstrj')
    assert trj.natoms == 1000
    assert [trj.has_attr(x) for x in 'tqxvf'] == [False, False, True, False, True]

    trj.read_frame()
    assert trj.f.shape == (1000, 3)
    assert abs(np.abs(trj.f).sum() - 28545.625) < 0.01

def test_reader(datafile):
    reader = TrajReader(datafile('methanol_1728_cg.trr'), skip=50, every=10, frames=25)
    nread = 0
    sumx = 0.0
    
    while reader.next_frame():
        nread += 1
        sumx += reader.traj.x.sum()
    
    assert nread == 25
    assert abs(sumx - 3195925.3359375) < 0.01
