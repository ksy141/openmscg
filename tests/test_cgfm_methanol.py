import pytest
from mscg.cli import cgfm



def test_1s(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("methanol_1728.data"),
        traj    = datafile("methanol_1728_cg.trr"),
        names   = 'MeOH',
        cut     = 8.0,
        pair    = ['model=BSpline,type=MeOH:MeOH,min=2.8,resolution=0.2'],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [25.981301124765647, 22.881032170827996, 13.88184871755871, 8.01630720872137, 2.115753165031487, 0.964077302114605, -1.662418639009231, 0.10251662595348932, 2.36532693124631, 0.8071321899945998, 0.782147198440079, 0.24986837054257383, 0.2124712298774945, -0.03528114161418795, -0.0221500717375667, -0.04352751810321006, -0.0013434035836443101, -0.07298427692488398, -0.062417399788417136, -0.07481467427180885]
    
    print(coeffs)
    
    for i in range(20):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff))
        assert abs(diff)<0.01



def test_2s(datafile):
    
    coeffs = cgfm.main(
        top     = datafile("methanol_1728_2s.data"),
        traj    = datafile("methanol_1728_2s.trr,frames=500"),
        names   = 'CH3,OH',
        cut     = 8.0,
        pair    = ['model=BSpline,type=CH3:CH3,min=2.9', 'model=BSpline,type=CH3:OH,min=2.8', 'model=BSpline,type=OH:OH,min=2.5'],
        bond    = ['model=BSpline,type=CH3:OH,min=1.35,max=1.65,resolution=0.01'],
        verbose = 0,
        save    = 'return'
    )
    
    benchmark = [20.2301991456678, 18.59956636643768, 17.98744867497155, 13.61571330246726, 11.98390192398072, 8.28587922329998, 6.51494860984118, 4.671167578717215, 3.6537890388041427, 2.650483368349131, 1.9134969346078208, 1.3780471469904896, 0.9041407183069098, 0.6299941709708665, 0.3607671058174531, 0.22596977136463717, 0.1411253321493879, 0.062309264443229435, 0.133945922796477, -0.05291872101137389]
    
    print("")
    
    for i in range(20):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff))
        assert abs(diff)<0.01
