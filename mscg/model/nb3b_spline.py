from mscg import *
import numpy as np
from ..core import cxx_model_nb3b_bspline as lib

class Nb3bSpline(Model):

    def __init__(self, **kwargs):
        self.min = -0.9
        self.max =  0.9
        self.order = 6
        self.resolution = 0.5
        self.gamma_ij = self.gamma_ik = 1.2
        self.a_ij = self.a_ik = 4.0
        
        self.serialized_names = ['min', 'max', 'resolution', 'order'] + \
            ['gamma_ij', 'a_ij', 'gamma_ik', 'a_ik']

        super().__init__(**kwargs)
        self._h = lib.create(self.min, self.max, self.resolution, self.order)

    def setup(self, top, pairlist):

        self.nparam = lib.get_npars(self.min, self.max, self.resolution, self.order)
        super().setup(top, pairlist)

        lib.setup(self._h, self.tid, pairlist._h, self.dF, self.dU)
        lib.setup_ex(self._h, self.tid_ij, self.gamma_ij, self.a_ij,
            self.tid_ik, self.gamma_ik, self.a_ik)

    def compute_fm(self):
        self.dF.fill(0)
        lib.compute_fm(self._h)

    def compute_rem(self):
    	assert 0

    def compute_table(self, x, force=True):
        vals = np.zeros(x.shape[0])
        lib.get_table(self._h, self.params, x, vals)
        return vals
