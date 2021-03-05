import pytest
from mscg.cli import cged

def test_rpc1(datafile):
    
    chi, groups = cged.main(
        pc    = datafile("edcg_rpc1_pc.npy"),
        ev    = datafile("edcg_rpc1_ev.npy"),
        npc   = 24,
        sites = 10,
        save  = 'return'
    )
    
    ref_chi = 179.92686
    ref_groups = [[  0,  44], [ 45,  89], [ 90, 135], [136, 192], [193, 238], 
                  [239, 285], [286, 295], [296, 321], [322, 358], [359, 371]]

    assert(abs(chi - ref_chi)<1.0)
    assert((groups == ref_groups).all())
    