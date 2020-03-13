#

import os
import numpy as np
from .core import cxx_traj as lib


class Trajectory:
    
    def __init__(self, filename, fmt=''):
        if not os.path.isfile(filename):
            raise Exception('File not found: ' + filename)
        
        if fmt == '':
            fmt = filename.split('.')[-1]
        
        fmt = fmt.lower()
        
        if fmt == 'trr':    
            self.h = lib.open_trr(filename)
        elif fmt == 'lammpstrj':
            self.h = lib.open_lmp(filename)
        else:
            raise Exception('Unknown file format: ' + fmt)
        
        st = lib.get_status(self.h)
        
        if st > 0:
            err_msg = [ "ready",
                "Failed to open the file",
                "Invalid header or format",
                "No frame data"
            ]
            raise Exception('Error: ' + err_msg[st])
        
        self.natoms = lib.get_natoms(self.h)
        self.has_force = lib.has_force(self.h)
        self.box = np.array(lib.get_box(self.h))
        
        self.x = None
        self.f = None
    
    def __del__(self):
        lib.close(self.h)
    
    def get_status(self):
        return lib.get_status(self.h)
    
    def rewind(self):
        lib.rewind(self.h)
    
    def read_frame(self):
        if lib.next_frame(self.h):
            self.box = np.array(lib.get_box(self.h))
            self.x = np.array(lib.get_x(self.h)).reshape(self.natoms, 3)

            if self.has_force:
                self.f = np.array(lib.get_f(self.h)).reshape(self.natoms, 3)
            
            return True
        else:
            return False
        
        
