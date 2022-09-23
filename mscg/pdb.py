from mscg import *
import pandas as pd
import numpy as np

class PDB:
    
    def __init__(self, rows):
        
        atom_list = []
        
        def _getint(row, start, end):
            s = row[start:end].strip()
            
            try:
                i = int(s)
                return i
            except:
                return np.nan
        
        for irecord, row in enumerate(rows):
            if len(row)<6:
                continue
            
            header = row[:6].strip()
            
            if header == 'ATOM' or header =='HETATM':
                atom_list.append({
                    'RECORD':  irecord + 1,
                    'SERIAL':  _getint(row, 6, 11),
                    'NAME':    row[12:16].strip(),
                    'RESNAME': row[17:20].strip(),
                    'CHAIN':   row[21].strip(),
                    'RESSEQ':  _getint(row,22,26),
                    'X':       float(row[30:38].strip()),
                    'Y':       float(row[38:46].strip()),
                    'Z':       float(row[46:54].strip()),
                    'ELEMENT': row[76:78].strip()
                })
            
            else:
                header = "UNKNOWN"
                
        self.atoms = pd.DataFrame(atom_list).set_index('RECORD')
        
    def get_coords(self):
        return self.atoms[['X', 'Y', 'Z']].values
    
    
    @staticmethod
    def read(filename):
        with open(filename, 'r') as f:
            return PDB(f.read().strip().split("\n"))