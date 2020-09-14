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
        
        types = self.type.split(',')
        
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
        
        self.dF = np.zeros(shape=(top.n_atom * 3, self.nparam))
        self.dU = np.zeros(self.nparam)

class ModelGroup:
    
    def __init__(self):
        self.items = []
        
    def empty(self):
        self.items = []
    
    def __iadd__(self, other):
        self.items.append(other)
        return self
    
    def compute_fm(self):
        for model in self.items:
            model.compute_fm()
    
    def compute_rem(self):
        for model in self.items:
            model.compute_rem()

models = ModelGroup()

class ModelArgAction(argparse.Action):
    
    def __init__(self, option_strings, dest, **kwargs):
        if dest not in Model.styles.keys():
            raise Exception('Illegal model style: [%s]' % (dest))
        
        super(ModelArgAction, self).__init__(option_strings, dest, **kwargs)
    
    def __call__(self, parser, namespace, values, option_string=None):
        if len(values)==1:
            values = values[0].split(' ')
        
        values = [v.split('=') for v in values]
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

"""
class ModelBase:
    
    STYLE_PAIR, STYLE_BOND, STYLE_ANGLE, STYLE_DIHEDRAL = 'Pair', 'Bond', 'Angle', 'Dihed'
    
    def __init__(self, types, top, geolist, style = STYLE_PAIR, nbody = 2, nparam = 2):
        if self.__class__.__name__ == 'ModelBase':
            raise TypeError("Cannot initiate the abstract base class for model.")
        
        self.types   = types
        self.top     = top
        self.geolist = geolist
        
        self.style   = style
        self.nbody   = nbody
        self.nparam  = nparam
        
        self.dudl = np.zeros(nparam)
        self.dudl2 = np.zeros(nparam)
        
        self.name = self.style + "_" + "_".join(self.types)
        
        if top is not None:
            get_type = getattr(top, "get_" + self.style.lower() + "_type")
            self.type_id = get_type(*(self.types))
        
        model.append(self)
    
    def check_top(self):
        if not isinstance(top, Topology):
            raise ValueError("Value passed to [top] is not an instance of class [Topology].")
    
    def check_list(self):
        if self.style == ModelBase.STYLE_PAIR:
            if not isinstance(geolist, PairList):
                raise ValueError("Value passed to [geolist] is not an instance of class [PairList].")
        else:
            if not isinstance(geolist, BondList):
                raise ValueError("Value passed to [geolist] is not an instance of class [BondList].")
    
    def require(self, func_name):        
        if getattr(self, func_name, None) is None:
            raise NotImplementedError("Required function [%s] is not implemented in this model." % (func_name))
    
    def set_params(self, params):
        self.params = params
        return self        
        
    def serialize(self, extra_data = {}):
        data = {
            'name':  self.name,
            'model': str(self.__class__.__name__),
            'style': self.style,
            'nbody': self.nbody,
            'types': self.types,
            'nparam': self.nparam,
            'dudl': self.dudl
        }
        
        data.update(extra_data)
        return data

class models:
    required_attrs = ['style', 'nbody', 'nparam', 'name']
    required_funcs = []
    
    items = []
    
    def __init__(self):
        raise TypeError("Cannot initiate an object from this class.")
        
    @classmethod
    def empty(cls):
        cls.items = []
    
    @classmethod
    def check(cls, item):
        if not issubclass(type(item), ModelBase):
            raise TypeError("Object is not a class inherited from the [ModelBase].")
        
        for k in cls.required_attrs + cls.required_funcs:
            if getattr(item, k, None) is None:
                raise NameError("Object is not from a fully implemented model class. Hint: missing [%s]" % (k))
    
    @classmethod
    def append(cls, item):
        cls.check(item)
        cls.items.append(item)
    
    @classmethod
    def call(cls, func):
        for m in cls.items: getattr(m, func)()
    
    @classmethod
    def get(cls, name):
        for m in cls.items:
            if m.name == name:
                return m
        
        return None
    
    @classmethod
    def extract(cls, attr):
        z = {}
        
        for m in cls.items:
            z[m.name] = getattr(m, attr)
        
        return z
        
    @classmethod
    def serialize(cls):
        z = {}
        
        for m in cls.items:
            z[m.name] = m.serialize()
        
        return z

import argparse
import importlib

class ModelCreator:
    def __init__(self, n, argname, args):
        segs = args.split(",")
        
        if len(segs) < n + 1:
            raise ValueError('incorrect number of fields for option --' + argname + ' ' + args)
        
        self.ntype   = n
        self.argname = argname
        self.style   = segs[0]
        self.types   = segs[1:n+1]
        self.kwargs  = {}
        
        for i in range(n+1, len(segs)):
            w = segs[i].split('=')
            self.kwargs[w[0]] = w[1]
    
    def create(self, top, geolist):
        model_module = importlib.import_module("mscg.model." + self.argname + "_" + self.style.lower())
        model_class = getattr(model_module, self.argname.capitalize() + self.style)
        model_class(self.types, top, geolist, **(self.kwargs))



def BuildModelArgAction(n, arg_name):
    class ModelArgAction(argparse.Action):
        nbody = n
        name = arg_name
        
        def __call__(self, parser, namespace, values, option_string=None):
            getattr(namespace, self.dest).append(ModelCreator(ModelArgAction.nbody, ModelArgAction.name, values))
            return
        
        def help():
            msg = "define new " + ModelArgAction.name + " model with format: "
            msg += ",".join(["style,type" + str(i+1) for i in range(ModelArgAction.nbody)]) 
            msg += ",kwargs"
            return msg
            
    return ModelArgAction
"""