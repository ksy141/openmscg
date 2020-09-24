"""
Class modules in this sub-package serve as the potential or force functions
to model the CG force-field. They can be used in all procedures developed
in this software, such as the **force-matching** and **relative-entroy** 
methods.

"""

from ..topology import Topology
from ..pairlist import PairList
from ..bondlist import BondList

import numpy as np
import argparse, importlib

class Model:
    styles = {
        'pair': 2,
        'bond': 2,
        'angle': 3,
        'dihedral': 4
    }
    
    def __init__(self, **kwargs):
        self.type = ""
        
        for k,v in kwargs.items():
            if not hasattr(self, k):
                raise KeyError('undefined argument name: ' + k)
            
            arg_type = type(getattr(self,k))
            
            try:
                setattr(self, k, arg_type(v))
            except:
                raise TypeError('unmatched type for argument: ' + k + ' ' + str(arg_type))
        
        self.nparam = 0
    
    def setup(self, top, itemlist):
        self.top = top
        self.list = itemlist
        
        types = self.type.split(':')
        
        if self.style == 'pair':
            self.tid = top.pair_tid(*types)
            self.name = '-'.join(types)
        else:
            if len(types)>1:
                self.tid = top.bonding_tid(self.style, types)
                self.name = '-'.join(types)
            else:
                self.tid = getattr(top, 'names_' + self.style).index(self.types[0])
                self.name = types[0]
        
        self.name = self.style.capitalize() + '_' + self.name
        self.dF = np.zeros(shape=(top.n_atom * 3, self.nparam))
        self.dU = np.zeros(self.nparam)
        self.params = np.zeros(self.nparam)

class ModelGroup:
    
    def __init__(self):
        self.items = []
        
    def empty(self):
        self.items = []
    
    def __iadd__(self, other):
        self.items.append(other)
        return self
    
    def __getitem__(self, key):
        if type(key)==int:
            return self.items[key]
        elif type(key)==str:
            for m in self.items:
                if m.name == key:
                    return m
        
        return None
    
    def compute_fm(self):
        for model in self.items:
            model.compute_fm()
    
    def compute_rem(self):
        for model in self.items:
            model.compute_rem()
    
    def serialize(self):
        serialized = {}
        
        for model in self.items:
            serialized[model.name] = {name:getattr(model, name) for name in ['style', 'type', 'params', 'nparam'] + model.serialized_names}
        
        return serialized
                

models = ModelGroup()

class ModelArgAction(argparse.Action):
    
    def __init__(self, option_strings, dest, **kwargs):
        if dest not in Model.styles.keys():
            raise Exception('Illegal model style: [%s]' % (dest))
        
        super(ModelArgAction, self).__init__(option_strings, dest, **kwargs)
    
    def __call__(self, parser, namespace, values, option_string=None):
        values = [v.split('=') for v in values.split(',')]
        
        if sum([1 if len(v)!=2 else 0 for v in values])>0:
            raise Exception('Illegal model arguments: [%s]' % (' '.join(values)))
        
        args = {v[0]:v[1] for v in values}
        
        model_name = args.pop('model', None);
        module_name = "mscg.model.%s_%s" % (self.dest, model_name.lower())
        model_module = importlib.import_module(module_name)
        
        model_class = getattr(model_module, self.dest.capitalize() + model_name)
        model = model_class(**args)
        setattr(model, 'style', self.dest)
        getattr(namespace, self.dest).append(model)
        
        global models
        models += model
        return
    
    @classmethod
    def help(cls, style):
        return "add a model declaration for %s-style interactions." % (style)
