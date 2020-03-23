#

from .core import cxx_bspline as lib

class BSpline:
    
    def __init__(self, order, resolution, xmin, xmax):
        self.h = lib.create(order, resolution, xmin, xmax)
        
    def __del__(self):
        lib.destroy(self.h)
    
    def interpolate(self, xmin, dx, n, coeffs):
        return lib.interp(self.h, xmin, dx, n, coeffs)