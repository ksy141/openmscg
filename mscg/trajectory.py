"""

A `trajectory` is a time-series data displaying changes of atomic 
properties from MD simulations. The most important and necessary 
atomic properties are the coordinates (*x*). Some trajectory data 
may contain other optional information such as force (*f*), 
velocities (*v*), atom types (*t*) or partial chargs (*q*).

Trajectories are processed by frame data, including the header and
body information. The header part usually contains the data about
the frame number, time and system box. The body part contains *Nx3*
matrices (for coordinates, velocities and forces) and/or *N*-length
arrays (for types and charges, etc).

Examples with Trajectory APIs
-----------------------------

Create a trajectory object for reading a Gromacs/TRR file::
    
    >>> from mscg import Trajectory
    >>> trj = Trajectory('tests/data/methanol_1728_cg.trr', fmt='trr')

Check number of atoms in this trjectory::
    
    >>> trj.natoms
    1728

Check available attributes in this trjectory::
    
    >>> [trj.has_attr(x) for x in 'tqxvf']
    [False, False, True, False, True]

Checkout the attribute arrays for a frame::
    
    >>> trj.read_frame()
    >>> trj.x
    [[40.93094  11.893124 40.36594 ]
     [39.87031  19.339375 19.890625]
     [46.330624 35.58281  29.184376]
     ...
     [21.187813 23.75219  40.038437]
     [24.632812  5.304374 35.20875 ]
     [39.371876 18.761564  4.821562]]
    >>> trj.v
    None
    
"""

import os
import numpy as np
from .core import cxx_traj as lib


class Trajectory:
    """
    A class providing low-level APIs to read or write trajectory data in files.
    
    Attributes:
        natoms : int
            Number of atoms
        box : numpy.array(shape=3, dtype=numpy.float32)
            a 3-element vector for dimensions of the simulation box
        t : numpy.array(shape=natoms, dtype=numpy.int32) or None
            an vector of atom type IDs, if avaialabe
        q : numpy.array(shape=natoms, dtype=numpy.float32) or None
            an vector of partial charges, if avaialabe
        x,v,f : numpy.array(shape=(natoms, 3), dtype=numpy.float32) or None
            a matrix of coordinates, velocities and forces, if avaialabe
    """
    
    _datadefs = {
        't': { 'matrix': False, 'dtype': np.int32 },
        'q': { 'matrix': False, 'dtype': np.float32 },
        'x': { 'matrix': True,  'dtype': np.float32 },
        'v': { 'matrix': True,  'dtype': np.float32 },
        'f': { 'matrix': True,  'dtype': np.float32 },
    }
    
    _data_names = 'tqxvf'
    
    def __init__(self, filename, mode='r', fmt=''):
        """Create a Trajectory object. In read-mode, an exception will be raised up if the file doesn't exist.
        
        :param filename: path and name of the trjectory file.
        :type filename: str
        
        :param mode: mode of the operation ("r" - read, "w" - new & write, or "a" - append), default to "r".
        :type mode: str
        
        :param fmt: format of the file, "trr" or "lammps".
        :type fmt: str
        
        :return: a `Trajectory` object
        :rtype: Trajectory()
        """
        
        self._h = None
               
        for name in type(self)._data_names:
            setattr(self, "_" + name, None)
        
        self._file = filename
        self._mode = mode
        self._data = {}
        
        if mode not in ["r", "w", "a"]:
            raise Exception('Unsupported mode: ' + mode)
                
        if fmt == '':
            fmt = filename.split('.')[-1]
        
        self.fmt = fmt.lower()
        
        if ("r" in mode) and (not os.path.isfile(filename)):
            raise Exception('File not found: ' + filename)
        
        if self.fmt == 'trr':    
            self._h = lib.open_trr(filename, mode)
        elif self.fmt == 'lammpstrj':
            self._h = lib.open_lmp(filename, mode)
        else:
            raise Exception('Unsupported file format: ' + self.fmt)
        
        self._status = lib.get_status(self._h)
        
        if self._status > 0:
            err_msg = [ "ready",
                "Failed to open the file",
                "Invalid header or format",
                "No frame data"
            ]
            raise Exception('Error: ' + err_msg[st])
        
        if "r" in mode:
            self.natoms = lib.get_natoms(self._h)
            self.box = np.ndarray(shape=(3), dtype=np.float32)
            self.__allocate__()            
    
    def __allocate__(self):
        for attr,defs in type(self)._datadefs.items():
            shape = (self.natoms, 3) if defs['matrix'] else (self.natoms,)
            self._data[attr] = np.zeros(shape=shape, dtype=defs['dtype']) if self.has_attr(attr) else None
        
    def __del__(self):
        if self._h is not None:
            lib.close(self._h)
    
    def __getattr__(self, name):
        if name in type(self)._data_names:
            return self._data[name]
        
        raise Exception('Unknown attribute: Trajectory.' + name)
        
    def __setattr__(self, name, value):
        if name in type(self)._data_names:
            self._data[name] = value
            return
        
        super().__setattr__(name, value)
        
    def has_attr(self, attr):
        attrs = lib.has_attr(self._h, attr)
        return attrs
    
    def rewind(self):
        """
        Rewind the file handle to the starting position.
        """
        lib.rewind(self._h)
    
    def read_frame(self) -> bool:
        """
        Read a frame from the trajectory file. Return *False* if file end is reached.
        
        :return: success to read a frame or not
        :rtype: bool
        """
        
        if self._mode != "r":
            raise Exception("File is not in reading mode.")
        
        if not lib.read_frame(self._h):
            return False
        
        if lib.get_natoms(self._h) != self.natoms:
            self.__allocate__()
                
        attrs = [self._data[attr] for attr in type(self)._data_names]
        lib.get_frame(self._h, self.box, *attrs)
        
        return True
    
    def write_frame(self):
        """
        Write a frame to the trajectory file. The current data stored in *t, q, x, v, f* will
        be dumped to the file, if the data is both avaiable in the object and supported by the 
        file format.
        
        :return: success to write a frame or not
        :rtype: bool
        """
        
        if self.mode not in ["w", "a"]:
            raise Exception("File is not in writing mode.")
        
        return lib.write_frame(self._h, self.box.astype(np.float32), 
            *[self._data[attr] for attr in type(self)._data_names])
            
    @classmethod
    def wrap_molecule(cls, x, box):
        """
        Wrap coordinates subject to the coordinates of the first row in *x*.
        This is a type of PBC correction, which updates the coordinates in
        *x[1:,:]* with distance shorter than half of the box dimension to
        *x[0,:]*.
        
        :param x: a group of coordinates to be wrapped. The first row is the reference.
        :type x: numpy.array(shape=(,3), dtype=numpy.float32)
        
        :param box: box dimensions
        :type box: numpy.array(shape=3, dtype=numpy.float32)
        
        :return: wrapped coordinates
        :rtype: numpy.array(shape=3, dtype=numpy.float32)
        """
        
        x_wrap = x - x[0]
        x_scaled = x_wrap/box
        x_wrap += (x_scaled>0.5) * box * -1.0
        x_wrap += (x_scaled<-0.5) * box
        return x_wrap + x[0]
    
    @classmethod
    def pbc(cls, box, x):
        """
        Wrap coordinates into the box.
        
        :param box: box dimensions
        :type box: numpy.array(shape=3, dtype=numpy.float32)
        
        :param x: a group of coordinates to be wrapped. The first row is the reference.
        :type x: numpy.array(shape=(,3), dtype=numpy.float32)
        
        :return: wrapped coordinates
        :rtype: numpy.array(shape=3, dtype=numpy.float32)
        """
        
        x = x/box
        x -= np.floor(x)
        x -= np.floor(x)
        return x * box
        






