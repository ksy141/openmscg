
from mscg import *

class ForceHarmonic:
    def __init__(self, params):
        self.params = params
    
    def setup(self, top):
        n = getattr(top, 'n_' + self.style)
        names = getattr(top, 'names_' + self.style)
        
        self.k = np.zeros(n)
        self.v0 = np.zeros(n)
        
        for i, t in enumerate(getattr(top, 'types_' + self.style)):
            name = names[t]
            self.k[i] = self.params[name][0]
            self.v0[i] = self.params[name][1]
    
    def compute(self, n, types, scalars, U, dU):
        dv = scalars - self.v0
        U[:] = 0.5 * self.k * dv * dv
        dU[:] = - self.k * dv
        

class Bond_Harmonic(ForceHarmonic):
    pass

class Angle_Harmonic(ForceHarmonic):
    pass

class Dihedral_Harmonic(ForceHarmonic):
    pass