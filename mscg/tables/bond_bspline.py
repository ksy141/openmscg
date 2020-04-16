
from . import tables
from ..core import cxx_table_bond_bspline as lib
import numpy as np

class TableBondBSpline:
    
    def __init__(self, blist, name, type_id, order=4, resolution=0.05, min=1.0, max=2.0):
        self.name = 'Bond_' + name
        self.h = lib.create(blist.h, type_id, order, resolution, min, max)
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
