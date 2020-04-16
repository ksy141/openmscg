
from . import tables
from ..core import cxx_table_angle_bspline as lib
import numpy as np

class TableAngleBSpline:
    
    def __init__(self, blist, name, type_id, order=4, resolution=3, min=0, max=180):
        self.name = 'Angle_' + name
        
        min = min / 180.0 * np.pi
        max = max / 180.0 * np.pi
        resolution = resolution / 180.0 * np.pi
        
        self.h =  lib.create(blist.h, type_id, order, resolution, min, max)
        tables.add(self)
    
    def __del__(self):
        lib.destroy(self.h)
        
    def setup_cache(self, ddx_factor=0.001):
        lib.setup_cache(self.h, ddx_factor)
    
    def compute(self):
        lib.compute(self.h)
    
    def get_spline(self):
        return lib.get_spline(self.h)
    
    def dump(self, xmin, dx, n):
        x = [xmin + dx*i for i in range(n)]
        f = lib.dump(self.h, xmin, dx, n)
        return np.array([x, f])
