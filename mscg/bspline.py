#

from .core import cxx_bspline as lib
import math

class BSpline:
    
    @classmethod
    def ncoeff(cls, order, resolution, xmin, xmax):
        return int(math.ceil((xmax - xmin)/resolution) + 1) + order - 2
    
    def __init__(self, order, resolution, xmin, xmax):
        self.h = lib.create(order, resolution, xmin, xmax)
        
    def __del__(self):
        lib.destroy(self.h)
    
    def interpolate(self, xmin, dx, n, coeffs):
        return lib.interp(self.h, xmin, dx, n, coeffs)