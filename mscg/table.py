#

import numpy as np

class LammpsTable:
    
    def __init__(self, model):
        self.model = model
    
    def dump(self, xmin, xmax, xinc):
        
        x = np.arange(xmin, xmax, xinc)
        u, f = self.model.eval_etbl(x)
        n = u.shape[0]
        
        txt = "# Table %s: id, r, potential, force\n\n" % (self.model.name)
        txt += "_".join(self.model.name.split("_")[1:]) + "\n"
        
        txt += "N %d R %f %f\n\n" % (n, x[0], x[-1])

        for i in range(n):
            txt += "%d %f %f %f\n" % (i+1, x[i], u[i], f[i])

        with open(self.model.name + ".table", "w") as f:
            f.write(txt)