import pytest
from mscg import *
from mscg.cli import cgfm

def test_bond(datafile, show):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        bond    = ["model=BSpline,type=1:1,min=2.8,max=3.2,resolution=0.05"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [38.10980703593296, 39.72307422887177, 34.04565534070925, 29.798003623715736, 21.971013893242038, 13.506968886269405, 4.345854769783856, -4.36361946932223, -13.450051593350974, -22.04370399181429, -29.65709458006238, -34.22932434725673, -39.44262519914467, -39.31793042646144]

    show("test_bond", coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01


def test_angle(datafile, show):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        angle   = ["model=BSpline,type=1:1:1,min=80,max=100,resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [0.5735653495532164, 0.5287794226283061, 0.5176558463503989, 0.31193537985451286, 0.1576141575238138, -0.12704771531551345, -0.27964323386430906, -0.45891761720568514, -0.5457420213332002, -0.6195491324447984]
    
    show("test_angle", coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01

def test_bond_angle(datafile, show):
    
    coeffs = cgfm.main(
        top     = datafile("4_site.top"),
        traj    = datafile("4_site_d_non_periodic.lammpstrj"),
        cut     = 15.0,
        bond    = ["model=BSpline,type=1:1,min=2.75,max=3.25,resolution=0.05"],
        angle   = ["model=BSpline,type=1:1:1,min=79,max=104,resolution=5"],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [62.18592913281134, 45.018147612208885, 43.650766154323456, 37.81880054371472, 30.439505878757824, 19.80019969520647, 10.195324118779503, -0.051290200176430024, -9.928157703913032, -20.00548225238321, -29.905710984148232, -38.33928612204292, -44.64822648237569, -42.41712596883157, -52.461729306266, 0.7074543255911069, 0.5863739691202232, 0.5366086600173553, 0.3624048642413299, 0.17642158446931008, -0.09364817889270535, -0.327872750594599, -0.5501688926593241, -0.7237024227054131, -0.6817899879321104, -0.9855831680005782]
    
    show("test_bond_angle", coeffs)
    print("")
    
    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01

def test_bond_angle_dihedral(datafile, show):
    
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
    
    benchmark = [61.721591722736775, 45.09864076926988, 43.577456062095735, 37.83112964212978, 30.406448624772395, 19.788450514423346, 10.179613892713462, -0.07085190819964282, -9.93974488895875, -20.03416305134067, -29.900095934018225, -38.37899854674515, -44.63478161522994, -42.47824568880378, -52.46695478121046, 0.8411191219407126, 0.7169876005836864, 0.6211513518900915, 0.46689490199499595, 0.24783797301016186, -0.028526836164756375, -0.30163910167544117, -0.5084357085794373, -0.7410894741231502, -0.656398456577513, -1.008836712717489, -5.264126921013661e-13, 7.942312215892399e-12, -2.3466332133737734e-11, -2.673002461418836e-12, 2.9818925456649736e-13, 1.2795473625527176e-12, 2.7737297279971926, -0.15033404955038066, 0.04457750269700256, -0.054471281901154356, -0.04372939577056112, -0.09267773837180471, -0.11366229866451683, -0.09525303937297736, -0.11594393824501026, -0.23523106512315728, -0.23506114968221148]
    
    show("test_bond_angle_dihedral", coeffs)
    print("")
    
    for i in range(len(benchmark)):
        if abs(benchmark[i])>0.1:
            diff = (coeffs[i] - benchmark[i]) / benchmark[i]
            print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
            assert abs(diff)<0.01
