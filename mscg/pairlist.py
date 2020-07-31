#

import numpy as np
from .core import cxx_pairlist as lib

class PageIterator:
    
    def __init__(self, h, page_size=50000):
        self.h = h
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
        
        self.inext, self.n = lib.fill_page(self.h, self.type_id, self.inext, self.page_size,
                                           self._index, self._vec, self._r)
        
        if self.n == 0:
            raise StopIteration
        
        self.index = None if self._index is None else self._index[:,:self.n]
        self.vec = None if self._vec is None else self._vec[:,:self.n]
        self.r = None if self._r is None else self._r[:self.n]        
        return self
    
    def __call__(self, type_id=0, index=False, vector=False, scalar=True):
        self.type_id = type_id
        self._index = self.index_cache if index else None
        self._vec = self.vec_cache if vector else None
        self._r = self.r_cache if scalar else None
        return self
    

class PairList:
    
    def __init__(self, top, page_size = 50000):
        self.h = lib.create(top.h, top.natoms)
        self.pages = PageIterator(self.h, page_size)
    
    def __del__(self):
        lib.destroy(self.h)
        
    def init(self, cut = 10.0, binsize = 5.0):
        self.cut = cut
        
        if binsize>cut:
            binsize = cut * 0.5
            
        lib.init(self.h, cut, binsize)
        
    def setup_bins(self, traj):
        lib.setup_bins(self.h, traj.h)
    
    def build(self, traj, reset_bins = False):
        return lib.build(self.h, traj.h, reset_bins)
    
    def get_pairs(self, start, count):
        return lib.get_pairs(self.h, start, count)
    
    def update_types(self, ntypes, types):
        return lib.update_types(self.h, ntypes, np.array(types, dtype=np.int32))
    
    
