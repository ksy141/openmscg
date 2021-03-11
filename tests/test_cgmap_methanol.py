import pytest
from mscg import Mapper
from mscg.cli import cgmap
import numpy as np


def test_mapper(datafile):
    
    box = np.array([1.0, 1.0, 1.0])
    x = np.array([[0.5, 0.5, 0.05], [0.6, 0.6, 0.95]])
    
    mapper = Mapper({
        'site-types': {
            'CG': {
                'index': [0, 1],
                'x-weight': [4.0, 6.0],
                'f-weight': [1.0, 1.0]
            }
        },
        'system': [{'anchor': 0, 'repeat': 1, 'offset': 0, 'sites': [['CG', 0]]}]
    })
    
    X, F = mapper.process(box, x, None)
    diff = X - np.array([0.56, 0.56, 0.99])
    assert abs(diff.dot(diff.T)[0])<0.0001

def test_methanol_2s(datafile):
    
    frames = cgmap.main(
        traj = datafile("methanol_1728_aa.trr,every=100"),
        map  = datafile("methanol_1728_2s_map.yaml"),
        out  = 'return'
    )
    
    ref = np.array([
        [36.40839869, 13.04069239, 40.63038457],
        [36.73720102, 16.65953782, 33.59948371],
        [47.13768972, 23.04829915, 23.09552498],
        [ 2.82534173, 21.65635748, 22.46876548],
        [ 3.2128164,  36.93597641, 29.73447329]
    ])
    
    for i in range(5):
        diff = frames['x'][i][i] - ref[i]
        assert abs(diff.dot(diff.T))<0.0001
