
from . import tables
from ..core import cxx_table_angle_bspline as lib
import numpy as np

class TableAngleBSpline:
    
    def __init__(self, blist, type_id, order=3, resolution=3, xmin=0, xmax=180):
        
        xmin = xmin / 180.0 * np.pi
        xmax = xmax / 180.0 * np.pi
        resolution = resolution / 180.0 * np.pi
        
        self.h =  lib.create(blist.h, type_id, order, resolution, xmin, xmax)
        tables.add(self)
    
    def __del__(self):
        lib.destroy(self.h)
        
    def setup_cache(self, ddx_factor=0.001):
        lib.setup_cache(self.h, ddx_factor)
    
    def compute(self):
        lib.compute(self.h)
    
    def dump(self, xmin, dx, n):
        x = [xmin + dx*i for i in range(n)]
        f = lib.dump(self.h, xmin, dx, n)
        return np.array([x, f])
