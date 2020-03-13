#

from .core import cxx_pairlist as lib

class PairList:
    
    def __init__(self, top, traj):
        self.h = lib.create(top.h, traj.h)
    
    def __del__(self):
        lib.destroy(self.h)
        
    def init(self, cut = 10.0, binsize = 5.0):
        lib.init(self.h, cut, binsize)
        
    def setup_bins(self):
        lib.setup_bins(self.h)
    
    def build(self, reset_bins = False):
        return lib.build(self.h, reset_bins)
    
    def get_pairs(self, start, count):
        return lib.get_pairs(self.h, start, count)
