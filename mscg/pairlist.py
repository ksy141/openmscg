#

from .core import cxx_pairlist as lib

class PairList:
    
    def __init__(self, top):
        self.h = lib.create(top.h, top.natoms)
    
    def __del__(self):
        lib.destroy(self.h)
        
    def init(self, cut = 10.0, binsize = 5.0):
        self.cut = cut
        
        if binsize>cut:
            binsize = cut * 0.5
            
        lib.init(self.h, cut, binsize)
        
    def setup_bins(self, traj):
        lib.setup_bins(self.h, traj.h)
    
    def build(self, traj, reset_bins = False):
        return lib.build(self.h, traj.h, reset_bins)
    
    def get_pairs(self, start, count):
        return lib.get_pairs(self.h, start, count)
