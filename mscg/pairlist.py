"""
`PairList` provides APIs for a high-performance implementation (by C++)
for the `Verlet-List <https://en.wikipedia.org/wiki/Verlet_list>`__ and
`Cell Linked-List <https://en.wikipedia.org/wiki/Cell_lists>`__ algorithms.
It is helpful for any processing/analysis needs the information of all 
non-bounded pairs of particles of the system.

Construction of pair-list is based on a topology and coodinates from a
trajectory. So, the `PairList` module is usually used with the `Topology`
and `Trajectory` modules together. However, the module is designed to be
flexible enough to accept customized data structure as arguments.

.. note::
    
    The `PairList` provides the list of atom index in pairs within a
    given cut-off distance, but also the vector (dx, dy, dz) components
    and scalar (dr) values of each pair. So, it is a powerful tool used
    for any modeling work related to pair-wise interactions.
    

Construct Pairs
---------------

The `PairList` object can be created independently::
    
    >>> from mscg import PairList
    >>> plist = PairList(cut = 10.0, binsize = 5.0)

Besides the cut-off distance, there's another optional parameter `binsize`
that defines the size of *cells* in the linked-list algorithm. This parameter
affects the computational performance very much, however, is hard to determine.
A common choice is half of the cut-off distance.

After creation, the `init()` function needs to be called to connect the object
to the topology information::
    
    >>> plist.init(types_atom, exclusion_map)

The parameter `types_atom` is a numpy list of atom types of the system, while 
`exclusion_map` is a 2-D NumPy array defines the excluded pairs. If *None* is
passed, no exclusion is considered for the pair-list construction.

Before building up the pair-list, the **cells** (bins) need to setup by calling
the `setup_bins()` function::
    
    >>> plist.setup_bins(box)

The parameter `box` is a 3-element NumPy array defining the dimensions of the
simulation box. After setting up the cells, a 2-d NumPy array with shape (N,3)
storing the coordinates can be passed to the function `build()` to construct
the pair-list, where *N* is the total number of atoms::
    
    >>> plist.build(x)

While processing a trajectory, this function will be called for every frame of
data. The function `setup_bins()` is required to be called if the size of box
is changed. If there's a change of the topology, i.e., number of atoms, the
`init()` function needs to be called again.

Extract Pairs
-------------

The data in a `PairList` object can be extracted via the `PageIterator` class,
which has an instance in the `PairList` object::

    >>> for page in pairlist.pages(top.types_atom, type_i = 0, type_j = 0, index = True, vector = True, scalar = True):
    ...     print(page.index)
    ...     print(page.vec)
    ...     print(page.r)

The parameter `page_size` passed in during the creaton of a `PairList` object
defines the page size (number of pairs extracted from each page). The iterator
is initiated by an array of type IDs for all atoms in the system, and two type
IDs for the certain type of pairs that are to be extracted.

There are three more *boolean* arguments that define what types of data to be
stored in the pages:

* If `index` is *True*, there will be the `page.index` attribute that is an array in shape of (Mx2), where M is the number pairs and each column contains 2 atom index of a pair.
* If `vector` is *True*, there will be the `page.vec` attribute that is an array in shape of (Mx3), which gives the *R(i) - R(j)* on x, y and z axis.
* If `sacalar` is *True*, there will be the `page.r` attribute that stores the distances of all pairs.

"""

import numpy as np
from .core import cxx_pairlist as lib

class PairList:
    """
    Pair-List constructor the using Verlet-List algorithm.
    
    Attribute:
        page : mscg.PageInterator
            An iterator of pages contaning the data of pairs
    """
    def __init__(self, cut = 10.0, binsize = 5.0, max_pairs = 2000000, page_size = 1000000):
        """
        Create an object in `PairList` class.
        
        :param cut: cut-off distance to construct the pair-list, unit in *Angstrom*.
        :type cut: float
        
        :param binsize: size of the cells in the licked-list algorithm, unit in *Angstrom*.
        :type binsize: float
        
        :param max_pairs: maximum number of pairs can be contained, default to 2,000,000.
        :type max_pairs: int
        
        :param page_size: number of pairs extracted to each page, default to 50,000.
        :type page_size: int
        """
        
        self.cut = cut
        
        if binsize>cut:
            binsize = cut * 0.5
        
        self._h = lib.create(cut, binsize, max_pairs)
        self._scalar_t, self._scalar_r = None, None
        self.pages = PageIterator(self._h, page_size)
    
    def __del__(self):
        lib.destroy(self._h)
        
    def init(self, types_atom, exmap = None):
        '''
        Initialize the memory to build up cell-lists. Needs to be called whenever the
        number of atoms or type-list of atoms is changed.
        
        :param types_atom: type-list of atoms
        :type types_atom: numpy.ndarray(N, dtype=np.int32)
        
        :param exmap: exclustion-map for pairs with N rows (N is the total number of atoms). Each row must be ended by -1.
        :type exmap: numpy.ndarray(shape=(N, M))
        '''
        
        self._types = types_atom.astype(np.int32).copy()
        
        if exmap is not None:
            if type(exmap) != np.ndarray or exmap.ndim!=2 or exmap.shape[0]!=types_atom.shape[0]:
                raise Exception('Illegal argument for [exmap].')
            
            self._exmap = exmap.copy()
        else:
            self._exmap = None
        
        lib.init(self._h, self._types, self._exmap)
        
    def setup_bins(self, box:np.ndarray):
        '''
        Setup the cells (bins) for the linked-list algorithm. Needs to be called whenever
        the PBC box is initiated or changed.
        
        :param box: dimensions of the system box (PBC).
        :type box: numpy.array(shape=(3,), dtype=float)
        '''
        lib.setup_bins(self._h, box)
        
    def build(self, x:np.ndarray):
        '''
        Build the pair-list.
        
        :param x: coordinates of the system.
        :type box: numpy.array(shape=(N,3), dtype=float)
        '''
        return lib.build(self._h, x)
    
    def update_types(self, t:np.ndarray):
        if t.dtype != np.int32:
            t = t.astype(np.int32)
        
        lib.update_types(self._h, t)
    
    def num_pairs(self):
        return lib.num_pairs(self._h)
    
    def get_scalar(self):
        n = self.num_pairs()
        
        if self._scalar_t is None or self._scalar_t.shape[0] < n:
            self._scalar_t = np.zeros(n, dtype=np.int32)
            self._scalar_r = np.zeros(n, dtype=np.float32)
        else:
            self._scalar_t.fill(0)
            self._scalar_r.fill(0)
            
        lib.get_scalar(self._h, self._scalar_t, self._scalar_r)
        return self._scalar_t, self._scalar_r

    @classmethod
    def get_tid(cls, i:int, j:int):
        return lib.get_tid(i, j)
    
class PageIterator:
    """
    Iterator class to extract data of pairs.
    
    Attributes:
        index : numpy.array(dtype=int) with shape=Mx2, or None
            atom index of pairs
        vec : numpy.array(dtype=float) with shape=Mx3, or None
            vectors (dx, dy, dz) for each pair
        r : numpy.array(dtype=float) with shape=M, or None
            distance of each pair of atoms
    """
    
    def __init__(self, h, page_size=1000000):
        self._h = h
        self.page_size = page_size
        
        self.index_cache = np.empty((2, page_size), dtype=int)
        self.vec_cache = np.empty((3, page_size))
        self.r_cache = np.empty(page_size)
        
    def __iter__(self):
        self.inext = 0
        self.n = self.page_size
        return self
    
    def __next__(self):
        if self.n < self.page_size:  
            raise StopIteration
        
        self.inext, self.n = lib.fill_page(
            self._h, self.type_id, self.inext, self.page_size,
            self._index, self._vec, self._r)
        
        if self.n == 0:
            raise StopIteration
        
        self.index = None if self._index is None else self._index[:,:self.n]
        self.vec = None if self._vec is None else self._vec[:,:self.n]
        self.r = None if self._r is None else self._r[:self.n]        
        return self
    
    def __call__(self, type_id, index=False, vector=False, scalar=True):
        """
        Initiate the iterator.
        
        :param type_id: type for the pairs of atoms to be extracted.
        :type type_id: int
        
        :param index: if to extract index data, default to be *False*.
        :type index: bool
        
        :param vector: if to extract vector data, default to be *False*.
        :type vector: bool
        
        :param scalar: if to extract scalar data, default to be *True*.
        :type scalar: bool
        
        """
        self.type_id = type_id
        self._index = self.index_cache if index else None
        self._vec = self.vec_cache if vector else None
        self._r = self.r_cache if scalar else None
        return self
