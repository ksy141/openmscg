#

from mscg import R2D
from ..table import *
import mscg.core.cxx_nb3b as lib

class NB3B_SW:
    def __init__(self, params):

        self.params = params
        self._h = None

    def __del__(self):
        if self._h is not None:
            lib.destroy(self._h)

    def setup(self, top):
        i = top.names_atom.index(self.params['type_i'])
        j = top.names_atom.index(self.params['type_j'])
        k = top.names_atom.index(self.params['type_k'])

        self._h = lib.create_sw(i, j, k,
            self.params['lambda'], self.params['cos0'],
            self.params['gamma_ij'], self.params['a_ij'],
            self.params['gamma_ik'], self.params['a_ik'])

    def compute(self, plist, U, dU, F):
        lib.compute_sw(self._h, plist._h, U, dU, F)
