class tables:
    
    all = []
    
    @staticmethod
    def add(tbl):
        tables.all.append(tbl)
    
    @staticmethod
    def compute_all():
        for table in tables.all:
            table.compute()
        
from .pair_bspline import *
from .bond_bspline import *
from .angle_bspline import *