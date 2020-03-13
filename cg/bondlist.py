#

from .core import cxx_bondlist as lib

class BondList:
    
    def __init__(self, top, traj):
        self.h = lib.create(top.h, traj.h)
    
    def __del__(self):
        lib.destroy(self.h)
    
    def build(self):
        return lib.build(self.h)
    
    def get_bonds(self):
        return lib.get_bonds(self.h)
    
    def get_angles(self):
        return lib.get_angles(self.h)
