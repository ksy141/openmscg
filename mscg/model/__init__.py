

from ..topology import Topology
from ..pairlist import PairList
from ..bondlist import BondList

import numpy as np

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

class model:
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