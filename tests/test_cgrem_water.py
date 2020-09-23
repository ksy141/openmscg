import pytest
from mscg.cli import cgderiv, cgrem

def test_water(datafile):

    cgderiv.main(
        top     = datafile("rem_water/system.data"),
        traj    = datafile("rem_water/reference.lammpstrj"),
        names   = 'SL',
        cut     = 7.0,
        pair    = ["model=BSpline type=SL,SL min=2.4 resolution=0.2 order=4"],
        save    = datafile("rem_water/model_ref"),
        verbose = 0
    )
    
    import os
    
    if os.path.exists(datafile("rem_water/restart.p")):
        os.remove(datafile("rem_water/restart.p"))
    
    args = """
--cut 7.0
--names SL
--pair model=BSpline type=SL,SL min=2.4 resolution=0.2 order=4
--save return
--verbose 0
"""
    args += "\n--top " + datafile("rem_water/system.data")
    args += "\n--traj " + datafile("rem_water/dump.lammpstrj")

    with open(datafile("rem_water/cgderiv.sh"), "w") as f:
        f.write(args)
    
    kwargs = {
        'ref'         : datafile("rem_water/model_ref.p"),
        'cgderiv-arg' : datafile("rem_water/cgderiv.sh"),
        'model'       : datafile("rem_water/model.txt"),
        'md'          : datafile("rem_water/md.inp"),
        'restart'     : datafile("rem_water/restart"),
        'table'       : datafile("rem_water/"),
        'optimizer'   : 'builtin,chi=0.5,t=298.15',
        'maxiter'     : 2,
        'verbose'     : 0
    }
    
    cgrem.main(**kwargs)
    
    import pickle
    result = pickle.load(open(datafile("rem_water/restart.p"), 'rb'))
    coeffs = result['iterations'][-1]['params']['Pair_SL-SL']
        
    benchmark = [4.929977,0.007959,-0.622645,-0.176770,0.111943,0.154674,0.120967,0.079169,0.041455,-0.001336,-0.025775,-0.036470,-0.024643,-0.001310,0.025274,0.057008,0.081885,0.079878,0.059818,0.038057,0.018955]
        
    print(','.join(["%0.6f" % (i) for i in coeffs]))
    print("")

    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01
    
    
    
