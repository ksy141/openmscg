#

from .core import cxx_matrix as lib

class Matrix:
    
    def __init__(self):
        self.h = lib.create()
        self.tables = []
        
    def __del__(self):
        lib.destroy(self.h)
        
    def add_tables(self, tables):
        for table in tables:
            self.tables.append(table)
            lib.add_table(self.h, table.h)
    
    def table_names(self):
        return [tbl.name for tbl in self.tables]
    
    def setup(self, natoms):
        lib.setup(self.h, natoms)
    
    def reset(self):
        lib.reset(self.h)
    
    def solve(self):
        lib.solve(self.h)
    
    def multiplyadd(self, traj):
        lib.multiplyadd(self.h, traj.h)
    
    def cov_X(self):
        return lib.cov_X(self.h)
    
    def cov_y(self):
        return lib.cov_y(self.h)
    
    def save(self, filename = "matrix"):
        X = self.cov_X()
        y = self.cov_y()
        
        tables = {}
        shift = 0
        
        for tbl in self.tables:
            param = tbl.get_spline()
            
            tables[tbl.name] = {
                'order':  param[0],
                'nbreak': param[1],
                'ncoeff': param[2],
                'xmin':   param[3],
                'xmax':   param[4],
                'dx':     param[5],
                'coeffs': y[shift:shift + param[2]]
            }
            
            shift += param[2]
        
        from datetime import datetime
        import pickle, os, sys, getpass, socket
        
        pickle.dump({
            'X':X, 'y':y, 'tables':tables,
            'os': sys.platform,
            'host': socket.gethostname(),
            'user': getpass.getuser(),
            'path': os.getcwd(),
            'time': str(datetime.now()),
            'file': filename + '.p'
        }, open(filename + '.p', 'wb'))
    
    
