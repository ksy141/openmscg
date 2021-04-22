import pytest
from mscg.cli import cgfm

def test_slab(datafile):
    kwargs = {
        'top'     : datafile("methanol_ucg.data"),
        'traj'    : datafile("methanol_slab_cg.trr") + ",every=5",
        'names'   : 'MeOH,Near,Far',
        'cut'     : 8.0,
        'pair'    : ['model=BSpline,type=Near:Near,min=3.0,resolution=0.05', 
                     'model=BSpline,type=Near:Far,min=3.0,resolution=0.05',
                     'model=BSpline,type=Far:Far,min=3.0,resolution=0.05'],
        'ucg'     : 'replica=25,seed=1234',
        'ucg-wf'  : 'RLE,I=MeOH,J=MeOH,high=Near,low=Far,rth=4.5,wth=3.5',
        'verbose' : 0,
        'save'    : 'return'
    }
    
    coeffs = cgfm.main(**kwargs)
    
    benchmark = [ 1.02706770e+01,  9.60111079e+00,  8.49875703e+00,  8.09731422e+00,
                  6.07569120e+00,  4.57822777e+00,  3.50977763e+00,  2.42697457e+00,
                  1.72940049e+00,  1.07372497e+00]
    
    print(coeffs)
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff))
        
        if abs(benchmark[i])>0.05:
            assert abs(diff)<0.01
        
        
