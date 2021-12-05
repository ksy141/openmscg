import numpy as np
import yaml
from mscg import Trajectory

class Mapper:
    
    def __init__(self, map_data):
        
        types = map_data['site-types']
        system = map_data['system']

        # build convert matrix

        for name, weights in types.items():
            v = np.array(weights['x-weight'])
            v /= v.sum()
            weights['x-weight'] = v[np.newaxis].T
            weights['f-weight'] = np.array(weights['f-weight'])[np.newaxis].T
    
        # unpack tree-like topology to flatten sites map

        def unpack_group(group, root_anchor=0):
            unit_sites = []
            
            if 'groups' in group:
                for item in group['groups']:
                    for site in unpack_group(item, root_anchor + group['anchor']):
                        unit_sites.append(site[:])
            
            if 'sites' in group:
                for site in group['sites']:
                    unit_sites.append(site[:])
            
            sites = []

            for i in range(group['repeat']):
                offset = root_anchor + group['anchor'] + i * group['offset']

                for site in unit_sites:
                    sites.append([site[0], site[1]+offset])

            return sites
        
        self.types = types
        self.sites = unpack_group({'anchor':0, 'offset':0, 'repeat':1, 'groups':system})
    
    def get_types(self):
        return [list(self.types.keys()).index(s[0]) + 1 for s in self.sites]
    
    def process(self, box, x, f):
        x_list = []
                
        for site_name, site_anchor in self.sites:
            site_index = [index + site_anchor for index in self.types[site_name]['index']]
                        
            x_aa = x[site_index].copy()
            x_aa = Trajectory.wrap_molecule(x_aa, box)
            
            x_cg = np.matmul(x_aa.T, self.types[site_name]['x-weight']).T
            x_list.append(x_cg)
        
        X = np.concatenate(x_list, axis=0)
        
        if f is not None:
            f_list = []
            
            for site_name, site_anchor in self.sites:
                site_index = [index + site_anchor for index in self.types[site_name]['index']]
                f_aa = f[site_index].copy()
                f_cg = np.matmul(f_aa.T, self.types[site_name]['f-weight']).T
                f_list.append(f_cg)
        
            F = np.concatenate(f_list, axis=0)
        else:
            F = None
        
        return Trajectory.pbc(box, X), F
    
    @classmethod
    def build_from_yaml(cls, yaml_file):
        with open(yaml_file, 'r') as f:
            map_data = yaml.load(f.read(), Loader=yaml.FullLoader)
        
        return Mapper(map_data)
