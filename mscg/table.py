#

import numpy as np

class Table:
    
    def __init__(self, model, force=True, prefix=''):
        self.model = model
        self.force = force
        self.prefix = ''
    
    def compute(self, xmin, xmax, xinc):
        self.x = np.arange(xmin, xmax, xinc)
        self.v = self.model.compute_table(self.x)
        self.n = self.v.shape[0]
        
        if(self.force):
            self.f = self.v
            u = np.zeros(self.n)
            u[:-1] = 0.5 * (self.f[:-1] + self.f[1:])
            u = np.cumsum(u[::-1])[::-1]
            self.u = u * xinc
            
        else:
            self.u = self.v
            f = np.zeros(self.n+1)
            f[1:-1] = -np.diff(self.u, 1) / xinc
            f[0] = 2.0 * f[1] - f[2]
            f[-1] = 2.0 * f[-2] - 2.0 * f[-3]
            self.f = 0.5 * (f[:-1] + f[1:])
                
    def dump_lammps(self, xmin, xmax, xinc):
        self.compute(xmin, xmax, xinc)
        
        txt = "# Table %s: id, r, potential, force\n\n" % (self.model.name)
        txt += self.model.name.split("_")[1] + "\n"
        txt += "N %d R %f %f\n\n" % (self.n, self.x[0], self.x[-1])
        
        for i in range(self.n):
            txt += "%d %f %f %f\n" % (i+1, self.x[i], self.u[i], self.f[i])

        with open(self.model.name + ".table", "w") as f:
            f.write(self.prefix + txt)