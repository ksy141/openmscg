#

from mscg import *
from .mm import *
from .core import cxx_force as lib

class Force:
    
    def __init__(self, funcs):
        self.funcs = []
        self.nfunc = {'pair':0, 'bond':0, 'angle':0, 'dihedral':0}
        
        for key, params in funcs.items():
            if key not in globals():
                raise Exception("Unknown force function: " + key)
            
            self.funcs.append(globals()[key](params))
            registered = False
            
            for style in self.nfunc:
                if key.startswith(style.capitalize()):
                    setattr(self.funcs[-1], 'style', style)
                    self.nfunc[style] += 1
                    registered = True
                    break
            
            if not registered:
                raise Exception("Invalid force funcation name: " + key)
        
    def setup(self, top, plist, blist):
        self.top = top
        self.plist = plist
        self.blist = blist
        
        for func in self.funcs:
            func.setup(top)
        
        self.f = np.zeros((top.n_atom, 3), dtype=np.float32)
        self.e = 0.0
        self.U = {k:None for k in self.nfunc.keys()}
        self.dU = {k:None for k in self.nfunc.keys()}
        screen.info("Force functions: " + str(self.nfunc))
        
    def compute(self):
        
        # Prepare U and dU arrays
        for style, nfunc in self.nfunc.items():
            if nfunc > 0:
                n = self.plist.num_pairs() if style=='pair' else self.blist.num_items(style)
            
                if self.U[style] is None or self.U[style].shape[0] < n:
                    self.U[style] = np.zeros(n, dtype=np.float32)
                    self.dU[style] = np.zeros(n, dtype=np.float32)
            
                self.U[style].fill(0.0)
                self.dU[style].fill(0.0)
        
        # calculate U and dU
        for func in self.funcs:
            if func.style == 'pair':
                types, R = self.plist.get_scalar()
                func.compute(self.plist.num_pairs(), types, R, self.U['pair'], self.dU['pair'])
            else:
                types = getattr(self.blist, func.style + '_type')
                S = self.blist.get_scalar(func.style)
                func.compute(self.blist.num_items(func.style), types, S, self.U[func.style], self.dU[func.style])
        
        # calculate F
        self.f.fill(0.0)
        
        if self.nfunc['pair']>0:
            lib.compute_pair(self.plist._h, self.dU['pair'], self.f)
        
        if self.nfunc['bond']>0:
            lib.compute_bond(self.blist._h, self.dU['bond'], self.f)
            
        if self.nfunc['angle']>0:
            lib.compute_angle(self.blist._h, self.dU['angle'], self.f)
        
        if self.nfunc['dihedral']>0:
            lib.compute_dihedral(self.blist._h, self.dU['dihedral'], self.f)
        
        # calculate E
        self.e = 0.0
        
        for style, nfunc in self.nfunc.items():
            if nfunc > 0:
                self.e += self.U[style].sum()
        
        return self.f




