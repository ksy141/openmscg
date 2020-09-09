import numpy as np
import time

class UCG:
    
    weighting_funcs = []
    replica = 100
    seed = int(time.time())
    
    top = None
    traj = None
    weights = None
    
    @classmethod
    def init(cls, plist, blist):
        for weighting in cls.weighting_funcs:
            setattr(weighting, "plist", plist)
            setattr(weighting, "blist", blist)
    
    @classmethod
    def process(cls):
        for weighting in cls.weighting_funcs:
            weighting.compute(cls.top, cls.traj, cls.weights)
        

import argparse
import importlib

class UCGArgAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        args = values.split(',')
        sf_name = args[0]
        kwargs = {}
        
        for arg in args:
            k, v = arg.split('=')
            
            if not hasattr(UCG, k):
                raise Exception("Wrong parameter name for UCG: " + k)
            
            try:
                vtype = type(getattr(UCG, k))
                setattr(UCG, k, vtype(v))
            except:
                raise Exception("Incorrect format for UCG parameter: " + arg)
        
        np.random.seed(UCG.seed)
        return

    def help():
        msg = "settings for UCG modeling"
        return msg

class WFArgAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        
        args = values.split(',')
        wf_name = args[0]
                                
        wf_module = importlib.import_module("mscg.ucg.wf_" + wf_name.lower())
        wf_class = getattr(wf_module, "Weighting" + wf_name)
        wf = wf_class()
        
        for arg in args[1:]:
            k, v = arg.split('=')
            
            if not hasattr(wf, k):
                raise Exception("Incorrect option for UCG weighting function %s: %s" % (wf_name, k))
            
            try:
                vtype = type(getattr(wf, k))
                setattr(wf, k, vtype(v))
            except:
                raise Exception("Incorrect format for UCG weighting function parameter: " + arg)
        
        wf.init()
        UCG.weighting_funcs.append(wf)
        return

    def help():
        msg = "define new state-function for UCG"
        return msg
    

class UCGSpawner:
    def __init__(self, top, traj):
        if UCG.weighting_funcs == []:
            UCG.replica = 1
        else:
            UCG.top = top
            UCG.traj = traj
            UCG.weights = [[] for _ in range(top.natoms)]
            UCG.process()
        
    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i >= UCG.replica:
            raise StopIteration
        else:
            self.i += 1
        
        if UCG.weighting_funcs == []:
            return None
        else:
            
            r = np.random.random(UCG.top.natoms)
            types = []
            
            for i in range(UCG.top.natoms):
                p = 0.0
                
                for state in UCG.weights[i]:
                    p += state[1]
                    if r[i]<=p:
                        types.append(state[0])
                        break
                
                if len(types)<=i:
                    types.append(UCG.weights[-1][0])
            
            return [UCG.top.get_atom_type(t) for t in types]



    