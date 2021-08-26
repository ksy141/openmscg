from mscg import *
import numpy as np

class PairYukawa(Model):
    
    def __init__(self, **kwargs):
        self.n = 1
        self.kappa = '1.0'
        self.qfile = 'none'
        self.q = None
        
        self.serialized_names = ['n', 'kappa', 'qfile']
        super().__init__(**kwargs)
        
    def setup(self, top, pairlist):
        if self.qfile != 'none':
            with open(self.qfile, "r") as f:
                self.q = np.array([float(_) for _ in f.read().strip().split("\n")])
        
        if self.q.shape[0] != top.n_atom:
            raise Exception("Incorrect number of charges in the qfile.")
                    
        self.nparam = self.n
        super().setup(top, pairlist)
        
        self.kappa = [float(x) for x in self.kappa.split(':')]
        
    def compute_rem(self):
        self.dU.fill(0)
        
        for page in self.list.pages(self.tid, index=True):
            qa = np.ones(page.r.size)
            qb = np.ones(page.r.size)
            
            if self.q is not None:
                for i in range(page.r.size):
                    qa[i] = self.q[page.index[0][i]]
                    qb[i] = self.q[page.index[1][i]]
            
            self.dU += np.array([np.sum(qa * qb * np.exp(-k * page.r) / page.r) for k in self.kappa])
                    
    def compute_table(self, x, force=True):
        if force:
            raise Exception("Cannot generate force table for the Gauss/Cut model.")
        
        vals = np.zeros(x.shape[0])
        
        for i, k in enumerate(self.kappa):
            vals += self.params[i] * np.exp(-k * x) / x
        
        return vals
        
