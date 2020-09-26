import pytest
from mscg.cli import cgderiv

def test_1frame(datafile):

    mean, var = cgderiv.main(
        top     = datafile("water_512_cg.data"),
        traj    = datafile("water_512_cg.lammpstrj"),
        names   = 'SL',
        cut     = 7.0,
        pair    = ["model=BSpline,type=SL:SL,min=2.4,max=7.0,resolution=0.2,order=4"],
        verbose = 0,
        save    = 'return'
    )

    benchmark = [0.841046,25.837122,186.406502,392.214687,287.086246,200.801884,205.038367,248.046392,292.732254,350.580963,386.895521,448.283933,489.094345,543.298264,539.376924,570.928150,582.165001,572.471255,656.737124,721.033952,795.402621,888.816173,955.052062,753.485567,531.901492,271.472152]
    
    coeffs = mean['Pair_SL-SL']
    print(','.join(["%0.6f" % (i) for i in coeffs]))
    print("")

    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01
