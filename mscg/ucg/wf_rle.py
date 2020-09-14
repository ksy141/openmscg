#!/user/bin/env python3

import numpy as np

class WeightingRLE:
    
    def __init__(self, **kwargs):
        self.target = "undefined"
        self.low    = "undefined"
        self.high   = "undefined"
        self.rth    = 0.0
        self.wth    = 0.0
        return
    
    def init(self):
        pass
    
    def compute(self, top, traj, weights):
        pair_type = top.pair_tid(self.target, self.target)
        rho = np.zeros(top.n_atom)
        
        for page in self.plist.pages(pair_type, index=True):
            w = 0.5 * (1.0 - np.tanh((page.r - self.rth) / (0.1 * self.wth)))
            
            for i in range(w.shape[0]):
                rho[page.index[0][i]] += w[i]
                rho[page.index[1][i]] += w[i]
        
        p = 0.5 * (1.0 + np.tanh((rho - self.wth) / (0.1 * self.wth)))
        
        for i in range(top.n_atom):
            weights[i] = [(self.high, p[i]), (self.low, 1.0-p[i])]
        