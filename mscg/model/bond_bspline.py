from mscg import *
import numpy as np
from ..core import cxx_model_bond_bspline as lib

class BondBSpline(Model):
    
    def __init__(self, **kwargs):
        self.min = 0.8
        self.max = 1.6
        self.order = 6
        self.resolution = 0.1
        super().__init__(**kwargs)
        
    def setup(self, top, bondlist):
        self.nparam = lib.get_npars(self.min, self.max, self.resolution, self.order)
        super().setup(top, bondlist)
        self._h = lib.create(self.min, self.max, self.resolution, self.order, self.tid, bondlist._h, self.dF, self.dU)
        lib.setup_cache(self._h, self.resolution * 0.01)
    
    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)
        
    def compute_rem(self):
        raise Error('not impletmented')