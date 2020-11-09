import pytest
from mscg.cli import cgfm

def test_bond(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        bond    = ["model=BSpline,type=1:1,min=2.8,max=3.2,resolution=0.05"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [39.85091642,37.97855976,34.87891868,29.44192467,22.12103692,13.44038147,4.38742087,-4.39991485,-13.40422791,-22.12920724,-29.48013058,-34.61095986,-38.68172789,-40.36834862]

    print(coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01


def test_angle(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        angle   = ["model=BSpline,type=1:1:1,min=80,max=100,resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [0.59065808, 0.5332647, 0.51386409, 0.31500864, 0.15369519, -0.12284413, -0.28469427, -0.45104363, -0.5593834, -0.58713964]
    
    print(coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01

def test_bond_angle(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        bond    = ["model=BSpline,type=1:1,min=2.75,max=3.25,resolution=0.05"],
        angle   = ["model=BSpline,type=1:1:1,min=79,max=104,resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [51.95082208296276, 50.746028061110565, 41.720861268243425, 38.39308071446055, 30.26278804868963, 19.860813277824406, 10.173976336889353, -0.036132367084263706, -9.931535050591954, -20.00493926124949, -29.888470122092315, -38.37831464880469, -44.49624483740308, -42.946177289001525, -51.50452873124365, 0.6634766557253255, 0.609532860111007, 0.526643592082018, 0.367746155922438, 0.172700862490772, -0.09024431519117204, -0.331366243035184, -0.544056478083025, -0.7332084784523047, -0.650113855350348, -1.077520884421709]
    
    print(", ".join([str(c) for c in coeffs]))
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01
