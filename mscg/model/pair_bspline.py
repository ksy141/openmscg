from mscg import *
import numpy as np
from ..core import cxx_model_pair_bspline as lib

class PairBSpline(Model):
    
    def __init__(self, **kwargs):
        self.min = 2.0
        self.order = 6
        self.resolution = 0.1
        super().__init__(**kwargs)
        
    def setup(self, top, pairlist):
        self.nparam = lib.get_npars(self.min, pairlist.cut, self.resolution, self.order)
        super().setup(top, pairlist)
        self._h = lib.create(self.min, pairlist.cut, self.resolution, self.order, self.tid, pairlist._h, self.dF, self.dU)
        lib.setup_cache(self._h, self.resolution * 0.01)        
    
    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)
        
    def compute_rem(self):
        self.dU.fill(0)
        lib.compute_rem(self._h)
        print(self.dU)
    
    @classmethod
    def compute_table(xmin, dx, n, params):
        pass
        