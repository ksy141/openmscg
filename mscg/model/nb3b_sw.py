from mscg import *
import numpy as np
from ..core import cxx_model_nb3b_sw as lib

class Nb3bSW(Model):

    def __init__(self, **kwargs):
        self.gamma_ij = self.gamma_ik = 1.2
        self.a_ij = self.a_ik = 4.0
        self.theta0 = 180.0
        
        self.serialized_names = ['gamma_ij', 'a_ij', 'gamma_ik', 'a_ik', 'theta0']

        super().__init__(**kwargs)
        self._h = lib.create(self.gamma_ij, self.a_ij, self.gamma_ik, self.a_ik, self.theta0)

    def setup(self, top, pairlist):
        self.nparam = 1
        super().setup(top, pairlist)
        lib.setup(self._h, self.tid, pairlist._h, self.dF, self.dU)
        lib.setup_ex(self._h, self.tid_ij, self.tid_ik)
   
    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)

    def compute_rem(self):
    	fatal("NB3B/SW cannot work for REM!")

    def compute_table(self, x, force=True):
        fatal("NB3B/SW is not a tabulated potential style!")
