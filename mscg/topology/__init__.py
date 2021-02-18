"""

A topology contains **three** groups of lists: `names` and `types` and `items`.
In each group, there are four lists, containing the information of
atoms and bonding structures: 2-body bonds, 3-body angles and
4-body dihedral torsions.

The `names` are used to define the names of types that are used to map
a set of parameters for calculating potential energies. For example,
the types of atoms decide the parameters of potential functions 
for pair-wise interations. 

The `types` are used to specify the type `id` for the bonding items in
the system. Each of the types list is a `numpy.array(dtype=int)`, in which
each element represents the index of the type in a types list.

The `items` contain the indices of atoms that form the bonds, angles ... in
the system. Each of the items list is a 2-D `numpy.ndarray(dtype=int)`, in which
each bond, angles ... are stored as columns.

The **Topology** class supports following features:
    
    * Manually building up topologies by adding type names, atoms and bonding items.
    * Replicating, combining topologies by using arithmetic operators `*` and `+`.
    * Reading topologies from standard MD input files.

Examples
--------

Create Topology from Scrach
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a new empty topology object::
    
    >>> from mscg import Topology
    >>> top = Topology
    
Add four names of atom types::
    
    >>> top.add_names('atom', ['CT','HT','OH','HO'])
    4

Add six atoms to build up a methanol molecule::
    
    >>> top.add_atoms(['CT'] + ['HT'] * 3 + ['OH', 'HO'])
    6
    >>> top.types_atom
    array([0, 1, 1, 1, 2, 3])

Get type names for atoms::
    
    >>> [top.names_atom[i] for i in top.types_atom]
    ['CT', 'HT', 'HT', 'HT', 'OH', 'HO']

Add bonds::
    
    >>> top.add_bondings('bond', ['?'] * 3, [[0, 0, 0,], [1,2,3]], True) # Add three C-H bonds
    3
    >>> top.add_bondings('bond', ['?'] * 3, [[0, 4], [4, 5]], True)
    5
    >>> top.ntype_bond
    3
    >>> top.types_bond
    array([0, 0, 0, 1, 2])
    >>> top.bond_atoms
    array([[0, 0, 0, 0, 4],
           [1, 2, 3, 4, 5]])
    >>> top.names_bond
    ['CT-HT', 'CT-OH', 'OH-HO']

Generate angles automatically::
    
    >>> top.gen_angles(True)
    >>> top.angle_atoms
    array([[1, 1, 1, 2, 2, 3, 5],
           [0, 0, 0, 0, 0, 0, 4],
           [2, 3, 4, 3, 4, 4, 0]])
    >>> top.names_angle
    ['HT-CT-HT', 'HT-CT-OH', 'CT-OH-HO']
    >>> top.types_angle
    array([0, 0, 1, 0, 1, 1, 2])
    >>> top.angle_atoms
    array([[1, 1, 1, 2, 2, 3, 5],
           [0, 0, 0, 0, 0, 0, 4],
           [2, 3, 4, 3, 4, 4, 0]])

Generate dihedral torsions automatically::

    >>> top.gen_dihedrals(True)
    >>> top.names_dihedral
    ['HT-CT-OH-HO']
    >>> top.types_dihedral
    array([0, 0, 0])
    >>> top.dihedral_atoms
    array([[1, 2, 3],
           [0, 0, 0],
           [4, 4, 4],
           [5, 5, 5]])

Check the linking map::
    
    >>> top.linking_map(True, False, False)
    array([[ 1,  2,  3,  4, -1],
           [ 0, -1, -1, -1, -1],
           [ 0, -1, -1, -1, -1],
           [ 0, -1, -1, -1, -1],
           [ 5,  0, -1, -1, -1],
           [ 4, -1, -1, -1, -1]])

Operations on Topologies
~~~~~~~~~~~~~~~~~~~~~~~~

Replicate content in a toplogy::
    
    >>> top.types_atom
    array([0, 1, 1, 1, 2, 3])
    >>> top *= 2
    >>> top.types_atom
    array([0, 1, 1, 1, 2, 3, 0, 1, 1, 1, 2, 3])
    
Combine topologies::
    
    >>> water = Topology('water')
    >>> water.add_names('atom', ['OW', 'HW'])
    2
    >>> water.add_atoms(['OW', 'HW', 'HW'])
    3
    >>> water.add_bondings('bond', ['?', '?'], [[0,0],[1,2]], True)
    2
    >>> water.gen_angles(True)
    >>> water.names_angle
    ['HW-OW-HW']
    >>> top = top + water * 3
    >>> top.names_angle
    ['HT-CT-HT', 'HT-CT-OH', 'CT-OH-HO', 'HW-OW-HW']
    >>> top.angle_atoms
    array([[ 1,  1,  1,  2,  2,  3,  ...,  13, 16, 19],
           [ 0,  0,  0,  0,  0,  0,  ...,  12, 15, 18],
           [ 2,  3,  4,  3,  4,  4,  ...,  14, 17, 20]])


Read Topologies from Files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Read the topology information from a LAMMPS data file::
    
    >>> top = Topology.read_file('tests/data/methanol_1728_2s.data')
    >>> top.reset_names(['CH3', 'OH'])
    >>> top.names_atom
    ['CH3', 'OH']
    >>> top.names_bond
    ['CH3-OH']
    >>> top.bond_atoms
    array([[   0,    2,    4, ..., 3450, 3452, 3454],
           [   1,    3,    5, ..., 3451, 3453, 3455]])

"""

import os, sys
from typing import Union
import numpy as np

from mscg.core import cxx_pairlist as core_pairlist

class Topology:
    """A class to represent the topology information of a molecular system.
    
    Names Attributes:
        names_atom : [str]
            Names of atom types.
        names_bond : [str]
            Names of bond types.
        names_angle : [str]
            Names of angle types.
        names_dihedral : [str]
            Names of dihedral types.
    
    Types Attributes:
        types_atom : numpy.array(dtype=int)
            Type IDs of atoms.
        types_bond : numpy.array(dtype=int)
            Type IDs of bonds.
        types_angle : numpy.array(dtype=int)
            Type IDs of angles.
        types_dihedral : numpy.array(dtype=int)
            Type IDs of dihedral torsions.
    
    Items Attributes:
        bond_atoms : numpy.ndarray(size=(2,n), dtype=int)
            Atom indices of bonds.
        angle_atoms : numpy.ndarray(size=(3,n), dtype=int)
            Atom indices of angles.
        dihedral_atoms : numpy.ndarray(size=(4,n), dtype=int)
            Atom indices of angles.
            
    .. note::
        
        All these attributes are copies of the topology information that are **NOT**
        for direct access. Modifications on these attributes directly won't change
        anything inside the `toplogy` object.
    
    The numbers of atoms, bonds, angles ... can be obtained from the lengths/shapes of
    the types attributes. However, this option will result extra work of making a copy.
    Following attributes are provided to get the counts more efficiently.
    
    Count Attributes:
        ntype_% : int
            number of type names for % (`%` should be 'atom', 'bond', 'angle' or 'dihedral')
        n_% : int
            number of atoms or bonding items for % (`%` should be 'atom', 'bond', 'angle' or 'dihedral')
            
    """
    
    bonded_types  = ['bond', 'angle', 'dihedral']
    bonded_natoms = {'bond':2, 'angle':3, 'dihedral':4}
    all_types     = ['atom'] + bonded_types
    
    def __setattr__(self, name, value):
        k = name.split('_')
        
        if k[0] in ['names', 'types', 'n', 'ntype']:
            raise Exception('Attempt to set a protected attribute: ' + name)
        
        super().__setattr__(name, value)
        
    def __getattr__(self, name):
        k = name.split('_')
        
        if k[0] in ['names', 'types']:
            _list = self.__dict__['_' + k[0]]
            
            if len(k)>1 and k[1] in _list:
                if k[0] == 'names':
                    return _list[k[1]].copy()
                else:
                    return np.array(_list[k[1]], dtype=np.int32)
        
        elif k[0] == 'n':
            if len(k)>1 and k[1] in type(self).all_types:
                return self._types[k[1]].shape[0]
        
        elif k[0] == 'ntype':
            if len(k)>1 and k[1] in type(self).all_types:
                return self._names[k[1]].__len__()
            
        elif len(k)>1 and k[1] == 'atoms':
            if k[0] in type(self).bonded_types:
                return self._bonded_atoms[k[0]].copy()
                    
        raise Exception('Unknown attribute: Topology.' + name)
        
    def __init__(self, system_name:str = 'undefined'):
        """
        Create an empty object in ``topology`` class. Initialize all type and item lists to be empty.
        
        :param system_name: name of the system, defaults to "undefined"
        :type system_name: str
        
        :return: an empty object in ``topology`` class.
        :rtype: Topology
        """
        
        self.system_name = system_name
        self._names = { t: [] for t in type(self).all_types }
        self._types = { t: np.array([], dtype=np.int32) for t in type(self).all_types }
        self._bonded_atoms = { t: np.ndarray((n, 0), dtype=np.int32) for t, n in type(self).bonded_natoms.items() }
            
    def copy(self):
        """
        Return a full (deep) copy of the current topology.
        """
        
        top = Topology()
        
        for t in type(self).all_types:
            top._names[t] = self._names[t][:]
            top._types[t] = self._types[t].copy()
        
        for t in type(self).bonded_types:
            top._bonded_atoms[t] = self._bonded_atoms[t].copy()
        
        return top
        
    def __mul__(self, n):
        """
        Repliacte the contents (atoms and bondings) for N times.
        """
        
        top = self.copy()
        shift = 0
                
        for i in range(1, n):
            top._types['atom'] = np.append(top._types['atom'], self._types['atom'])
            shift += self.n_atom
            
            for t in type(self).bonded_types:
                if self._bonded_atoms[t].shape[1] == 0:
                    continue
                    
                top._types[t] = np.append(top._types[t], self._types[t])
                
                atoms = self._bonded_atoms[t].copy()
                atoms += shift
                top._bonded_atoms[t] = np.append(top._bonded_atoms[t], atoms, axis=1)
        
        return top
    
    def __iadd__(self, o):
        """
        Merged by another topology object
        """
        
        assert type(o) == Topology
        
        for t in type(self).all_types:
            self.add_names(t, o._names[t])
        
        n = self.n_atom
        self.add_atoms([o._names['atom'][i] for i in o._types['atom']])
        
        for t in type(self).bonded_types:
            atoms = o._bonded_atoms[t].copy() + n
            atoms = [list(row) for row in list(atoms)]
            self.add_bondings(t, ['?'] * len(atoms[0]), atoms)
        
        return self
    
    def __add__(self, o):
        """
        Merge two topology objects into one
        """
        
        assert type(o) == Topology
        
        top = self.copy()
        top += o        
        return top
        
    def add_names(self, name_type:str, values) -> int:
        """
        Add a list of names to a type list.
        
        :param name_type: "atom", "bond", "angle" or "dihedral"
        :type system_name: str
        
        :param values: a list of names
        :type values: [str]
        
        :return: current number of type names in this list.
        :rtype: int
        
        Example::
            
            >>> top.add_names('atom', ['OW', 'HW])
            2
            >>> top.add_names('bond', ['OW-HW'])
            1
            
        """
        
        if name_type not in type(self).all_types:
            raise KeyError('Unknown name for a type list: ' + name_type)
        
        if type(values)!=str and type(values)!=list:
            raise ValueError('Name values must be in types of str or [str].')
        
        if type(values) == list:
            if sum([type(i)!=str for i in values]) > 0:
                raise ValueError('Name value is not set a str: ' + str([i for i in values if type(i)!=str]))
        else:
            values = [values]
        
        values = [t for t in values if t not in self._names[name_type]]
        self._names[name_type].extend(values)
        return self._names[name_type].__len__()
            
    def add_atoms(self, types:list) -> int:
        """
        Add a list of atoms defined by names of atom types. The names must 
        exist in the names list of atom types. If a type name cannot be found, 
        a `ValueError` will be raised up.
        
        :param types: a list of atom types
        :type types: [str]
                
        :return: current number of atom in the system.
        :rtype: int
        
        Example::
            
            >>> top.add_atoms(['OW', 'HW, 'HW])
            3
            
        """
        
        self._add_items_with_types('atom', types)
        return self._types['atom'].shape[0]
            
    def add_bondings(self, bonded_type:str, types:list, atoms:list, autotype=False) -> int:
        """
        Add a list of bonding items defined by types and index atoms. The type names must 
        exist in the current list. if the type is "**?**", the function will auto-generate the
        type name based on the types of bonded atoms. If `autotype` is **True**, unexisted
        type names will be automatically added to the name list. If `autotype` is **False** and 
        a type name cannot be found, a `ValueError` will be raised up.
        
        :param bonded_type: "bond", "angle" or "dihedral"
        :type: bonded_type: str
        
        :param types: a list of bond/angle/dihedral types
        :type types: [str]
        
        :param items: (MxN) 2-D list of atom indices, each column is the atoms for a bond, angle or dihderal.
        :type types: [[int], [int] ...]
        
        :param autotype: add unexisted type names automatically, default to *False*.
        :type autotype: bool
        
        :return: current number of items in the list for `bonded_type`.
        :rtype: int
        
        Example::
            
            >>> top.add_bondings('bond', ['?']*2, [[0, 0], [1, 2]], autotype=True)
            2
            
        """
        
        # check input data
        if bonded_type not in type(self).bonded_types:
            raise KeyError('Unknown name for a bonded type list: ' + bonded_type)
        
        nrow, ncol = len(atoms), len(atoms[0])
        
        if nrow != type(self).bonded_natoms[bonded_type]:
            raise ValueError('Incorrect count of rows for %s index atoms: %d(input) != %d(required)' % (bonded_type, nrow, type(self).bonded_natoms[bonded_type]))
            
        if sum([len(atoms[i])!=ncol for i in range(1, nrow)]) > 0:
            raise ValueError('Inconsistent counts of columns for %s index atoms.' % (bonded_type))
        
        atoms = np.array(atoms, dtype=np.int32)
        
        if (atoms>=len(self.types_atom)).sum() > 0 or (atoms<0).sum() > 0:
            raise ValueError('Incorrect index values for %s atoms.' % (bonded_type))
        
        # check types
        types = [self.bonding_name([self._types['atom'][atoms[j][i]]
                    for j in range(type(self).bonded_natoms[bonded_type])])
                        if types[i]=='?' else types[i] for i in range(ncol)]
        
        if autotype:
            for t in types:
                if t not in self._names[bonded_type]:
                    self._names[bonded_type].append(t)
                        
        # add index atoms
        self._bonded_atoms[bonded_type] = np.append(self._bonded_atoms[bonded_type], atoms, axis=1)
        
        # add types and return
        self._add_items_with_types(bonded_type, types)
        return self._types[bonded_type].__len__()
    
    def mask_bonds(self, mask:list) -> int:
        """
        Refine a new group of bonds from a subset of the existing ones by a mask.
        Angle and dihedrals are not affected. Needs to manually call gen_angels()
        and gen_dihedrals() if necessary.
        
        :param mask: list of boolean values.
        :type: mask: [boolean]
        
        :return: current number of items in the list for `bonded_type`.
        :rtype: int
        """
        
        if self._types['bond'].__len__() != len(mask):
            raise ValueError('Incorrect length of the mask.')
        
        self._types['bond'] = self._types['bond'][mask].copy()
        self._bonded_atoms['bond'] = self._bonded_atoms['bond'][:,mask].copy()
        return self._types['bond'].__len__()
    
    def bonding_name(self, atom_types) -> str:
        """
        Get the name for a group of bonding types.
        
        :param atom_types: type names of atoms in the bonding item.
        :type autotype: [str]
        
        :return: name of the bonding type.
        :rtype: str
        """
        
        if atom_types[0] > atom_types[-1]:
            atom_types = reversed(atom_types)
        
        return '-'.join([self._names['atom'][i] for i in atom_types])
    
    def pair_name(self, atom_types) -> (int, str):
        """
        Get the type id and name for a given pairs of atom types.
        
        :param atom_types: type IDs of atoms of the pair.
        :type autotype: [int]
        
        :return: ID and name of the pair type.
        :rtype: (int, str)
        """
        
        if atom_types[0] > atom_types[1]:
            atom_types = reversed(atom_types)
        
        return '-'.join([self._names['atom'][i] for i in atom_types])
    
    def _add_items_with_types(self, type_name:str, types):
        """
        Add a list of items defined by types.
        """
        
        try:
            index = [self._names[type_name].index(t) for t in types]
        except ValueError as e:
            raise Exception('Unknown %s-type name: %s' % (type_name, str(e)))
        
        self._types[type_name] = np.append(self._types[type_name], index)
    
    def linking_map(self, bond:bool = True, angle:bool = True, dihedral:bool = True) -> np.ndarray:
        """
        Get a linking map for atoms in the system. The linking map is returned as 2-D NumPy array in
        the shape (n,m), where *n* is the number of atoms and m is the maximum linked neighbors for
        atoms. For the row *i*, *map[i][j]* is the atom index of the *j-th* neighbor to atom *i*, while
        the empty cells are filled with -1.
        
        The function accepts three boolean arguments indicating the consideration of 1-2 atoms in bonds,
        1-3 atoms in angles or 1-4 atoms in dihedrals. The function can be used to construct the bonding
        map, as well as exclusion maps for pair-wise interactions.        
        
        :param bond: whether to exclude 1-2 bonded pairs, default to *True*.
        :type bond: bool
        
        :param angle: whether to exclude 1-3 pairs of angles, default to *True*.
        :type angle: bool
        
        :param dihedral: whether to exclude 1-3 pairs of dihedral tortions, default to *True*.
        :type angle: bool
        
        :return: 
        :rtype: numpy.array(shape=(natom,maxlink+1))
        """
        
        ex = [[] for _ in range(self.n_atom)]
        
        for t in type(self).bonded_types:
            if locals()[t]:    
                for ia in range(self.n_atom):
                    index = self._bonded_atoms[t]
                    n = type(self).bonded_natoms[t] - 1
                    ex[ia].extend(list(index[n][index[0] == ia]))
                    ex[ia].extend(list(index[0][index[n] == ia]))
        
        maxcols = max([len(i) for i in ex]) + 1        
        return np.array([i + [-1] * (maxcols - len(i)) for i in ex], dtype=np.int32)
        
    def reset_names(self, names):
        """
        Reset names for atom types and auto-generate names for bonding types.
        
        :param names: list of new names for atom types.
        :type names: [str]
        """
        
        if len(names) != self.ntype_atom:
            raise ValueException("inconsistent number of types for atoms.")
        
        m = {}
        
        for i in range(self.ntype_atom):
            m[self.names_atom[i]] = names[i]
        
        self._names['atom'] = names
        
        t = [one.split('-') for one in self.names_bond]
        self._names['bond'] = ['-'.join([m[tt[0]], m[tt[1]]]) for tt in t]
        
        t = [one.split('-') for one in self.names_angle]
        self._names['angle'] = ['-'.join([m[tt[0]], m[tt[1]], m[tt[2]]]) for tt in t]
        
        t = [one.split('-') for one in self.names_dihedral]
        self._names['dihedral'] = ['-'.join([m[tt[0]], m[tt[1]], m[tt[2]], m[tt[3]]]) for tt in t]
        
    def gen_angles(self, autotype=False):
        """
        Generate a new list of angles from existed bonds and replace the current one.
        
        :param autotype: add unexisted type names automatically, default to *False*.
        :type autotype: bool
        """
        
        links = self.linking_map(bond=True, angle=False, dihedral=False)
        angles = [[], [], []]

        for i in range(self.n_atom):
            link = list(links[i])
            n = link.index(-1) if -1 in link else len(link)

            for j in range(n-1):
                for k in range(j+1, n):
                    angles[0].append(link[j])
                    angles[1].append(i)
                    angles[2].append(link[k])
        
        self._types['angle'] = np.array([], dtype=np.int32)
        self._bonded_atoms['angle'] = np.ndarray((3, 0), dtype=np.int32)
        self.add_bondings('angle', ['?'] * len(angles[0]), angles, autotype)
        
    def gen_dihedrals(self, autotype=False):
        """
        Generate a new list of dihderals from existed bonds and replace the current one.
        
        :param autotype: add unexisted type names automatically, default to *False*.
        :type autotype: bool
        """
        
        links = self.linking_map(bond=True, angle=False, dihedral=False)
        diheds = [[], [], [], []]
        
        for ib in range(self.n_bond):
            ia, ja = self._bonded_atoms['bond'][0][ib], self._bonded_atoms['bond'][1][ib]
            ilink, jlink = list(links[ia]), list(links[ja])
            ni = ilink.index(-1) if -1 in ilink else len(ilink)
            nj = jlink.index(-1) if -1 in jlink else len(jlink)
            
            for k in range(ni):
                for l in range(nj):
                    if ilink[k] != ja and jlink[l] != ia:
                        diheds[0].append(ilink[k])
                        diheds[1].append(ia)
                        diheds[2].append(ja)
                        diheds[3].append(jlink[l])
        
        self._types['dihedral'] = np.array([], dtype=np.int32)
        self._bonded_atoms['dihedral'] = np.ndarray((4, 0), dtype=np.int32)
        self.add_bondings('dihedral', ['?'] * len(diheds[0]), diheds, autotype)
    
    def pair_tid(self, type_i:str, type_j:str) -> int:
        '''
        Get the internal type ID for a pair of atom types.
        
        :param type_i: name of type I
        :type type_i: str
        
        :param type_j: name of type J
        :type type_j: str
        
        :return: type ID of the pair.
        :rtype: int
        '''
        names = self._names['atom']
        return core_pairlist.get_tid(names.index(type_i), names.index(type_j))
    
    def bonding_tid(self, item, atom_names) -> int:
        '''
        
        '''
        
        name1 = '-'.join(atom_names)
        name2 = '-'.join(reversed(atom_names))
        
        if name1 in self._names[item]:
            return self._names[item].index(name1)
        elif name2 in self._names[item]:
            return self._names[item].index(name2)
        
        raise Exception('Unknown bonding type: ' + str(atom_names) + ' in ' + item)
    
    @classmethod
    def search_formatters(cls):
        formatters = {}
        import importlib
        
        def search_formatters(path):
            return [f[10:-3] for f in os.listdir(path) if f.startswith('formatter_') and f.endswith('.py')]
                
        for mod in search_formatters(os.path.dirname(__file__)):
            formatters[mod] = {'module': importlib.import_module("mscg.topology.formatter_" + mod) }
        
        for mod in search_formatters('.'):
            formatters[mod] = {'module': importlib.import_module('formatter_' + mod) }
                
        for name in formatters:
            if hasattr(formatters[name]['module'], 'inspect'):
                formatters[name]['inspector'] = getattr(formatters[name]['module'], 'inspect')
            
            if hasattr(formatters[name]['module'], 'read'):
                formatters[name]['reader'] = getattr(formatters[name]['module'], 'read')
            
            if hasattr(formatters[name]['module'], 'write'):
                formatters[name]['writer'] = getattr(formatters[name]['module'], 'write')
        
        return formatters
    
    def save(self, formatter, **kwargs):
        formatters = Topology.search_formatters()
        
        if formatter is None or formatter not in formatters:
            raise Exception("Unknown topology format.")
        
        if 'writer' not in formatters[formatter]:
            raise Exception("Cannot find a writer for the format: " + formatter)
        
        formatters[formatter]['writer'](self, **kwargs)
        
            
    @classmethod
    def read_file(cls, file:str, formatter=None):
        '''
        Read the topology from a file. If `formatter` is set to *None*, the format of
        file will be determined automatically.
        
        :param file: full path-name of the file
        :type file: str
        
        :param formatter: format of the file, *'lammps'* or *'cgtop'*.
        :type formatter: str
        
        :return: an object in ``topology`` class.
        :rtype: Topology
        '''
        if not os.path.isfile(file):
            raise Exception('File is not found: ' + file)
        
        with open(file, "r") as f:
            rows = [row.strip().split() for row in f.read().split("\n")]
        
        formatters = cls.search_formatters()
                        
        if formatter is None:            
            for name in formatters:
                if 'inspector' in formatters[name]:
                    if formatters[name]['inspector'](rows):
                        formatter = name
                        break
        
        if formatter is None:
            raise Exception("Unknown topology format.")
        
        if 'reader' not in formatters[formatter]:
            raise Exception("Cannot find a reader for the format: " + formatter)
        
        return formatters[formatter]['reader'](rows)
    
   
        
