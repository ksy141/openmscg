#

from .core import cxx_matrix as lib

class Matrix:
    
    def __init__(self):
        self.h = lib.create()
    
    def __del__(self):
        lib.destroy(self.h)
        
    def add_tables(self, tables):
        for table in tables:
            lib.add_table(self.h, table.h)
    
    def setup(self, natoms):
        lib.setup(self.h, natoms)
    
    def reset(self):
        lib.reset(self.h)
    
    def solve(self):
        lib.solve(self.h)
    
    def multiplyadd(self, traj):
        lib.multiplyadd(self.h, traj.h)
    
    
