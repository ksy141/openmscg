import pytest
from mscg import *

def test_methanol2s(datafile):
    top = Topology.read_file(datafile('methanol_1728_2s.data'))
    trj = Trajectory(datafile('methanol_1728_2s.trr'), fmt='trr')
    trj.read_frame()

    pairlist = PairList(page_size=500000)
    pairlist.init(top.types_atom, None)
    pairlist.setup_bins(trj.box)
    assert pairlist.build(trj.x) == 207591
    
    exmap = top.linking_map(True, True, True)
    pairlist.init(top.types_atom, exmap)
    pairlist.setup_bins(trj.box)
    assert pairlist.build(trj.x) == 205863

    for page in pairlist.pages(0, scalar=True, vector=True):
        assert abs(page.r.mean() - 7.592287814753148) < 0.001