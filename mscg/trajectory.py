#

import os
import numpy as np
from .core import cxx_traj as lib


class Trajectory:
    
    def __init__(self, filename, fmt=''):

        self.x = None
        self.f = None
        self.h = None
        
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
    
    def __del__(self):
        if self.h is not None:
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
        

class TrajReader:
    
    def __init__(self, file, skip, every, frames):
        
        segs = file.split(".")
        suffix = segs[-1] if len(segs)>1 else ""
        
        self.traj   = Trajectory(file, suffix)
        self.file   = file
        self.skip   = skip
        self.every  = every
        self.frames = frames
        self.nread  = 0
        
        for i in range(self.skip):
            self.traj.read_frame()
    
    def next_frame(self):
        
        if self.frames>0 and self.nread>=self.frames:
            return False
        
        for i in range(self.every-1):
            self.traj.read_frame()
        
        if self.traj.read_frame():
            self.nread += 1
            return True
        else:
            return False






