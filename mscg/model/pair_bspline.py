
from ..core import cxx_model_pair_bspline as lib
from . import ModelBase
from ..bspline import BSpline as bs
import numpy as np

class PairBSpline(ModelBase):
    
    def __init__(self, types, top, pairlist, **kwargs):
        
        self.order      = int(kwargs.get('order', 6))
        self.resolution = float(kwargs.get('resolution', 0.1))
        self.cache      = bool(kwargs.get('cache', False))
        self.min        = float(kwargs.get('min', 2.0))
        self.max        = pairlist.cut
        
        super().__init__(types, top, pairlist, style = ModelBase.STYLE_PAIR, nbody = 2,
                         nparam = bs.ncoeff(self.order, self.resolution, self.min, self.max))
        
        self.h = lib.create(pairlist.h, self.type_id, self.order, self.resolution, self.min)
        self.setup_cache()
    
    def __del__(self):
        return
        lib.destroy(self.h)
    
    def serialize(self):
        return super().serialize({
            'order'     : self.order,
            'resolution': self.resolution,
            'min'       : self.min,
            'max'       : self.max,
            'cache'     : self.cache
        })
    
    def eval_etbl(self, x:np.array, params = None):
        if params is not None:
            self.params = params
        
        xlist = list(x)
        xlist = [xlist[0] * 2.0 - xlist[1]] + xlist + [xlist[-1] * 2.0 - xlist[-2]]
        x = np.array(xlist)
        e = np.zeros(x.shape[0])
        lib.compute_etable(self.h, np.array(self.params), x, e)
        
        dx = np.diff(x)
        de = np.diff(e)
        f = - de / dx
        f = 0.5 * (f[:-1] + f[1:])
        return e[1:-1], f
    
    def eval_deriv_u(self):
        self.dudl = np.array(lib.compute_dudl(self.h))
    
    def setup_cache(self, ddx_factor=0.001):
        lib.setup_cache(self.h, ddx_factor)
        