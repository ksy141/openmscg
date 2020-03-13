
from . import tables
from ..core import cxx_table_bond_bspline as lib
import numpy as np

class TableBondBSpline:
    
    def __init__(self, blist, type_id, order=3, resolution=0.05, xmin=1.0, xmax=2.0):
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
