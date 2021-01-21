from mscg import *
import numpy as np

class PairGaussCut(Model):
    
    def __init__(self, **kwargs):
        self.n = 1
        self.rmh = '10.0'
        self.sigma = '1.0'
        self.serialized_names = ['n', 'rmh', 'sigma']
        super().__init__(**kwargs)
        
    def setup(self, top, pairlist):
        class GaussTerm:
            def __init__(self, rmh, sigma):
                self.rmh = rmh,
                self.factor1 = -1.0 / 2.0 / (sigma**2)
                self.factor2 = 1.0 / sigma / np.sqrt(np.pi * 2.0)
            
        self.nparam = self.n
        super().setup(top, pairlist)
        
        rmhs = [float(x) for x in self.rmh.split(':')]
        sigmas = [float(x) for x in self.sigma.split(':')]
        self.terms = [GaussTerm(rmhs[i], sigmas[i]) for i in range(self.n)]
        
    def compute_rem(self):
        self.dU.fill(0)
        
        for page in self.list.pages(self.tid):
            self.dU += np.array([np.sum(t.factor2 * np.exp(np.square(page.r - t.rmh) * t.factor1)) for t in self.terms])
                    
    def compute_table(self, x, force=True):
        if force:
            raise Exception("Cannot generate force table for the Gauss/Cut model.")
        
        vals = np.zeros(x.shape[0])
        
        for t in self.terms:
            vals += t.rmh * t.factor2 * np.exp(np.square(x - t.rmh) * t.factor1)
        
        return vals
        
