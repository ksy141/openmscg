#!/user/bin/env python3

import numpy as np

class WeightingRLE:
    
    def __init__(self, **kwargs):
        self.I      = "undefined"
        self.J      = "undefined"
        self.low    = "undefined"
        self.high   = "undefined"
        self.rth    = 0.0
        self.wth    = 0.0
        return
    
    def init(self):
        pass
    
    def compute(self, top, traj, weights):
        pair_type = top.pair_tid(self.I, self.J)
        tid = top.names_atom.index(self.I)
        types = top.types_atom
        rho = np.zeros(top.n_atom)
        
        for page in self.plist.pages(pair_type, index=True):
            w = 0.5 * (1.0 - np.tanh((page.r - self.rth) / (0.1 * self.rth)))
            
            for i in range(w.shape[0]):
                ia, ib = page.index[0][i], page.index[1][i]
                
                if types[ia] == tid:
                    rho[ia] += w[i]
                
                if types[ib] == tid:
                    rho[ib] += w[i]
        
        p = 0.5 * (1.0 + np.tanh((rho - self.wth) / (0.1 * self.wth)))        
        
        for i in range(top.n_atom):
            if types[i] == tid:
                weights[i] = [(self.high, p[i]), (self.low, 1.0-p[i])]
        
