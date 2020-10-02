from mscg import *
import numpy as np
from ..core import cxx_model_angle_bspline as lib

class BondBSpline(Model):
    
    def __init__(self, **kwargs):
        self.min = 0.0
        self.max = np.pi
        self.order = 6
        self.resolution = 0.1
        super().__init__(**kwargs)
        self._h = lib.create(self.min, self.max, self.resolution, self.order)
        lib.setup_cache(self._h, self.resolution * 0.01)
        
    def setup(self, top, bondlist):
        self.nparam = lib.get_npars(self.min, self.max, self.resolution, self.order)
        super().setup(top, bondlist)
        self.setup(self._h, self.tid, bondlist._h, self.dF, self.dU)
    
    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)
        
    def compute_rem(self):
        self.dU.fill(0)
        lib.compute_rem(self._h)
    
    def compute_table(self, x, force=True):
        vals = np.zeros(x.shape[0])
        lib.get_table(self._h, self.params, x, vals)
        return vals