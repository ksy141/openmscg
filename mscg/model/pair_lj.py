from . import ModelBase
import numpy as np

class PairLJ(ModelBase):
    
    def __init__(self, types, top = None, pairlist = None, **kwargs):        
        super().__init__(types, top, pairlist, style = ModelBase.STYLE_PAIR, nbody = 2, nparam = 2)
        
        self.alpha = int(kwargs.get('order_rep', 12))
        self.beta  = int(kwargs.get('order_att', 6))
    
    def __del__(self):
        pass
    
    def serialize(self):
        return super().serialize({
            'order_rep': self.alpha,
            'order_att': self.beta
        })
    
    def eval_etbl(self, x, params = None):
        if params is not None:
            self.params = params
        
        A = self.params[0]
        B = self.params[1]
        
        part1 = np.power(x, -self.alpha)
        part2 = np.power(x, -self.beta)
        e = A * part1 - B * part2
        
        part1 /= x
        part2 /= x
        f = A * self.alpha * part1 - B * self.beta * part2
        
        return e, f
    
    def eval_deriv_u(self):
        self.dudl.fill(0.0)
        
        for page in self.geolist.pages(self.type_id):
            self.dudl += np.array([np.power(page.r, -self.alpha).sum(), -np.power(page.r, -self.beta).sum()])
        
        
    
        
        