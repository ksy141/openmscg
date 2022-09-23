"""
Class modules in this sub-package serve as the potential or force functions
to model the CG force-field. They can be used in all procedures developed
in this software, such as the **force-matching** and **relative-entroy** 
methods. A model supported in MSCG usually provide three interfaces:

1. Calculate the force loadings on the parameters, which is used in FM.
2. Calculate the energy derivatives to the parameters, which is used in REM.
3. Generate tables by a given set of parameters.

"""

from ..topology import Topology
from ..pairlist import PairList
from ..bondlist import BondList

import numpy as np
import argparse, importlib
import os

class Model:
    """
    Base class for models.
    
    Attributes:
        style: str
            name of the style, `pair`, `bond`, `angle`, or `dihedral`
        type : str
            name of the targeted interaction type
        tid : int
            type ID of the targeted interaction defined in the topology referenced by `self.top`
        top : mscg.Topology
            associated topology object for this model
        list : mscg.PairList or mscg.BondList
            listing object for this model. If this model is for the `pair` style, this will be an object of `mscg.PairList`, otherwise `mscg.BondList`
        nparam: int
            number of parameters for this model
        params: numpy.array
            1-D array of parameters
        dU: numpy.array
            1-D array to store energy derivatives on parameters, with the shape of (1, self.nparam)
        dF: numpy.array
            2-D array to store force loadings of parameters on each dimension of atoms, with the shape of (self.top.n_atom x 3, self.nparam)
        serialized_names: [str]
            names of member attributes that are needed to stored during serialization
    """
    
    styles = {
        'pair': 2,
        'bond': 2,
        'angle': 3,
        'dihedral': 4,
        'nb3b': 3
    }
    
    def __init__(self, **kwargs):
        """
        Constructor of the model class. This function is used to parse the input
        argument ``**kwargs`` by extracting the values to set class attributes.
        
        The function will check the name and corresponding type is matched with
        existing member attribute in this class. So, this function is usually called
        at the end of the constructor in a child class, after declaring and setting
        default values of the attributes in the child class.
        """
        
        self.type = ""
        
        for k,v in kwargs.items():
            if not hasattr(self, k):
                raise KeyError('undefined argument name: ' + k)
            
            arg_type = type(getattr(self,k))
            
            try:
                setattr(self, k, arg_type(v))
            except:
                raise TypeError('unmatched type for argument: ' + k + ' ' + str(arg_type))
        
        self.module_name = type(self).__module__
        self.class_name = type(self).__name__
        self.nparam = 0
    
    def setup(self, top:Topology, itemlist):
        """
        Setup the model class, after receving the topology and listing objects. The model objects are constructed when parsing the runtime options, when the topology and listings may not be ready yet. Therefore, there need to be another function to setup the object after the topology and listing is created.
        
        One of the important role of this function is to allocate memory spaces (NumPy arrays) after receiving the number of atoms from the topology.        
        """
        
        self.top = top
        self.list = itemlist
        
        types = self.type.split(':')
        
        if self.style == 'pair':
            self.tid = top.pair_tid(*types)
            self.name = '-'.join(types)
        elif self.style == 'nb3b':
            self.tid = top._names['atom'].index(types[1])
            self.tid_ij = top.pair_tid(types[1], types[0])
            self.tid_ik = top.pair_tid(types[1], types[2])
            self.name = '-'.join(types)
        else:
            if len(types)>1: # this part is for ambigious A-B and B-A bonding types
                self.tid = top.bonding_tid(self.style, types)
                self.name = '-'.join(types)
            else: # for other styles, angle and dihedral
                self.tid = getattr(top, 'names_' + self.style).index(self.types[0])
                self.name = types[0]
        
        self.name = self.style.capitalize() + '_' + self.name
        self.dF = np.zeros(shape=(top.n_atom * 3, self.nparam))
        self.dU = np.zeros(self.nparam)
        self.params = np.zeros(self.nparam)
    
    def compute_fm(self):
        """
        Compute force loadings for the Force-Matching method, and stores the result in the attribute `self.dF`.
        """
        raise Error('Abstract method is not instantiated yet.')
    
    def compute_rem(self):
        """
        Compute energy derivatives for the Relative-Entropy method, and stores the result in the attribute `self.dU`.
        """
        raise Error('Abstract method is not instantiated yet.')
    
    def compute_table(self, x, force=True) -> np.array:
        """
        Compute tabulated values by given variable values.
        
        :param x: variable values for computing the table
        :type x: numpy.array
        
        :return: energy or force values
        :rtype: numpy.array
        """
        raise Error('Abstract method is not instantiated yet.')


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
            serialized[model.name] = {name:getattr(model, name) for name in ['module_name', 'class_name', 'style', 'type', 'params', 'nparam'] + model.serialized_names}
            
        return serialized

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


models = ModelGroup()
