import pytest
from mscg.cli import cgfm

def test_bond(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        bond    = ["model=BSpline type=1,1 min=2.8 max=3.2 resolution=0.05"],
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

'''
def test_angle(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        angle   = ["model=BSpline type=1,1,1 min=80 max=100 resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [0.6033390997989966, 0.5425583878460605, 0.3859922844497738, 0.1386004570285362, -0.10970664721238348, -0.339993126801044, -0.5268646865879502, -0.583997298766073]
    
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
        bond    = ["model=BSpline type=1,1 min=2.75 max=3.25 resolution=0.05"],
        angle   = ["model=BSpline type=1,1,1 min=79 max=104 resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [49.947571180460216, 46.57210091106926, 39.54644052444978, 30.325632914443283, 19.919855607160105, 10.114742565610829, 0.006338842257520952, -9.960716011681885, -19.977724701152578, -29.928635480479418, -40.79061530818936, -43.76065295897743, -51.48860803678937, 0.6651218666501006, 0.567848816733788, 0.37448211834472905, 0.06522166052678713, -0.23563017657869229, -0.5428967383628124, -0.760938588359891, -0.7783228220824304]
    
    print(coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01
'''
