from mscg import *
import numpy as np
from ..core import cxx_model_angle_bspline as lib

class AngleBSpline(Model):
    
    def __init__(self, **kwargs):
        self.min = 0.0
        self.max = 180.0
        self.order = 6
        self.resolution = 5.0
        self.serialized_names = ['min', 'max', 'resolution', 'order']
        super().__init__(**kwargs)
        self._h = lib.create(self.min * D2R, self.max * D2R, self.resolution * D2R, self.order)
        
        if bs_caching:
            lib.setup_cache(self._h, 0.001)
        
    def setup(self, top, bondlist):
        self.nparam = lib.get_npars(self.min * D2R, self.max * D2R, self.resolution * D2R, self.order)
        super().setup(top, bondlist)
        lib.setup(self._h, self.tid, bondlist._h, self.dF, self.dU)
    
    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)
        self.dF *= R2D
        
    def compute_rem(self):
        self.dU.fill(0)
        lib.compute_rem(self._h)
    
    def compute_table(self, x, force=True):
        vals = np.zeros(x.shape[0])
        lib.get_table(self._h, self.params, x * D2R, vals)
        return vals
