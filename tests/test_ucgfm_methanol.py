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
    
    benchmark = [ 9.76936203e+00,  9.67226099e+00,  8.10697754e+00,  8.09133619e+00,
                  5.83382416e+00,  4.48503194e+00,  3.31803244e+00,  2.31035101e+00,
                  1.58778322e+00,  9.49680476e-01,  6.05811162e-01,  1.00390383e-02]
    
    print(coeffs)
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff))
        
        if abs(benchmark[i])>0.05:
            assert abs(diff)<0.01
        
        
