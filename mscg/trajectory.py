#

import os
import numpy as np
from .core import cxx_traj as lib


class Trajectory:
    
    def __init__(self, filename, mode='r', fmt=''):

        self.x = None
        self.f = None
        self.h = None
        self.t = None
        
        self.file = filename
        self.mode = mode
        
        if mode not in ["r", "w", "a"]:
            raise Exception('Unsupported mode: ' + mode)
                
        if fmt == '':
            fmt = filename.split('.')[-1]
        
        self.fmt = fmt.lower()
        
        if ("r" in mode) and (not os.path.isfile(filename)):
            raise Exception('File not found: ' + filename)
        
        if self.fmt == 'trr':    
            self.h = lib.open_trr(filename, mode)
        elif self.fmt == 'lammpstrj':
            self.h = lib.open_lmp(filename, mode)
        else:
            raise Exception('Unknown file format: ' + self.fmt)
        
        st = lib.get_status(self.h)
        
        if st > 0:
            err_msg = [ "ready",
                "Failed to open the file",
                "Invalid header or format",
                "No frame data"
            ]
            raise Exception('Error: ' + err_msg[st])
        
        if "r" in mode:
            self.natoms = lib.get_natoms(self.h)
            self.box = np.ndarray(shape=(3), dtype=np.float32)
            self.x = np.ndarray(shape=(self.natoms,3), dtype=np.float32)
            self.t = np.ndarray(shape=(self.natoms,1), dtype=np.int32)
            self.v = np.ndarray(shape=(self.natoms,3), dtype=np.float32)
            self.f = np.ndarray(shape=(self.natoms,3), dtype=np.float32)
        
        
        
    def __del__(self):
        if self.h is not None:
            lib.close(self.h)
    
    def get_status(self):
        return lib.get_status(self.h)
    
    def has_type(self):
        return lib.has_force(self.h)
    
    def has_vel(self):
        return lib.has_vel(self.h)
    
    def has_force(self):
        return lib.has_force(self.h)
    
    def rewind(self):
        lib.rewind(self.h)
    
    def read_frame(self):
        if self.mode != "r":
            raise Exception("File is not in reading mode.")
        
        return lib.read_frame(self.h, self.box, self.t, self.x, self.v, self.f)
    
    def write_frame(self):
        if self.mode not in ["w", "a"]:
            raise Exception("File is not in writing mode.")
        
        return lib.write_frame(self.h, 
                               self.box.astype(np.float32), 
                               None if self.t is None else self.t.astype(np.int32), 
                               self.x.astype(np.float32), 
                               None if self.v is None else self.v.astype(np.float32), 
                               None if self.f is None else self.f.astype(np.float32))
            
    @classmethod
    def wrap_molecule(cls, x, box):
        x_wrap = x - x[0]
        x_scaled = x_wrap/box
        x_wrap += (x_scaled>0.5) * box * -1.0
        x_wrap += (x_scaled<-0.5) * box
        return x_wrap + x[0]
    
    @classmethod
    def pbc(cls, box, x):
        x = x/box
        x -= np.floor(x)
        x -= np.floor(x)
        return x * box
        






