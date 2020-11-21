import pytest
from mscg import *
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

def test_bond_angle_dihedral(datafile):
    
    coeffs = cgfm.main(
        top      = datafile("4_site.top"),
        traj     = datafile("4_site_d_non_periodic.lammpstrj"),
        cut      = 15.0,
        bond     = ["model=BSpline,type=1:1,min=2.75,max=3.25,resolution=0.05"],
        angle    = ["model=BSpline,type=1:1:1,min=77,max=104,resolution=5"],
        dihedral = ["model=BSpline,type=1:1:1:1,min=-60,max=60,resolution=10"],
        verbose  = 0,
        save     = 'return'
    )
    
    benchmark = [51.887739687974545, 50.67200515801255, 41.6712700971093, 38.363255446563926, 30.21912646069478, 19.822104887155774, 10.13230579500069, -0.07433130559691392, -9.973813862335781, -20.041303368228725, -29.936346085342052, -38.40576645774502, -44.5286241806391, -43.01615322146654, -51.59868161355428, 0.7559540517792762, 0.7417116309727021, 0.6078144129599137, 0.4682860495116117, 0.24117609402525275, -0.029224732482777682, -0.3073841929286729, -0.5073015712350326, -0.7507186112662794, -0.6346512533332529, -1.088630858922663, 0.17613322250449914, 0.16740138039262265, 0.24254503735823807, 0.08810822146061525, 0.12628926343945412, 0.09490495455348591, 0.05142623398802115, 0.037940422228853105, -0.0011573779623816338, -0.03386795501000489, -0.05889865138766803, -0.07846624836145266, -0.1358763247557846, -0.0678086695248048, -0.14986822165344638, -0.21193575005975673, -0.27729723183232324]
    
    print(", ".join([str(c) for c in coeffs]))
    print("")
    
    for i in range(len(benchmark)):
        if abs(benchmark[i])>0.1:
            diff = (coeffs[i] - benchmark[i]) / benchmark[i]
            print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
            assert abs(diff)<0.01
