"""
Similar to the `PairList` module, the `BondList` module provide a tool
calculate and store the geometric values for bonded interactions, including
bonds, angles and dihedral torsins.

`BondList` doesn't need to generate the list of bonding items, instead, feed
the items from a given topology object.

Usage
-----

The `PairList` object can be created independently::
    
    >>> from mscg import BondList
    >>> blist = BondList(bonds, angles, dihedrals)

Calculate values of bonding items::
    
    >>> blist.build(box, x)

Checkout scalar values of the bonds (angles or dihedral torsions)::
    
    >>> blist.get_scalar('bond')
    [1.5140172  1.1416421  1.1169515  1.4656354  1.2160345  0.92482877
     1.4008483  1.0002607  0.99983937 0.99994326 1.0006534  0.9999119
     0.9999411  0.999307   0.99996626 1.0001911  1.0003035  0.99956334
     0.9999446  1.0004575  0.99972695 1.0504069  1.4279559  1.0955634
     1.5082035  1.1287462  1.0932475 ]

"""

from .core import cxx_bondlist as lib
import numpy as np

class BondList:
        
    def __init__(self, bond_types, bonds, angle_types, angles, dihedral_types, dihedrals):
        """
        Create an object in BondList class.
        
        :param bonds: list of bonds (atom index) stored in columns
        :type bonds: numpy.array(shape=(2, N), dtype=numpy.int32)
        
        :param angles: list of angles (atom index) stored in columns
        :type angles: numpy.array(shape=(3, N), dtype=numpy.int32)
        
        :param dihedrals: list of dihderal torsions (atom index) stored in columns
        :type dihedrals: numpy.array(shape=(4, N), dtype=numpy.int32)
        """
        self._data = {}
        self._data['bond'] = bonds.transpose().copy()
        self._data['bond_type'] = bond_types.astype(np.int32)
        self._data['angle'] = angles.transpose().copy()
        self._data['angle_type'] = angle_types.astype(np.int32)
        self._data['dihedral'] = dihedrals.transpose().copy()
        self._data['dihedral_type'] = dihedral_types.astype(np.int32)
        
        self._h = lib.create(
            self._data['bond_type'], self._data['bond'], 
            self._data['angle_type'], self._data['angle'], 
            self._data['dihedral_type'], self._data['dihedral'])
    
    def __del__(self):
        if self._h is not None:
            lib.destroy(self._h)
        
    def build(self, box, x):
        """
        Calculate scalar and vector values for given bonds, angles and dihedral torsions.
        
        :param box: dimensions of the system PBC box 
        :type box: numpy.array(3, dtype=numpy.float32)
        
        :param x: coordinates of atoms
        :type x: numpy.array(shape=(N,3), dtype=numpy.float32)
        """
        return lib.build(self._h, box, x)
    
    def get_scalar(self, name):
        """
        Get scalar values from the list.
        
        :param name: *'bond'*, *'angle'* or *'dihedral'* 
        :type name: str
        
        :rtype: numpy.array(N, dtype=numpy.float32)
        """
        targets = ['bond', 'angle', 'dihedral']
        
        if name in targets:
            data = np.zeros(self._data[name].shape[0], dtype=np.float32)
            lib.get_scalar(self._h, targets.index(name), data)
            return data
