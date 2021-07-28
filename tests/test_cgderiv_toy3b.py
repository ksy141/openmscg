import pytest
from mscg.cli import cgderiv

def test_3site(datafile):

    mean, var = cgderiv.main(
        top     = datafile("toy_3site.top"),
        traj    = datafile("toy_3site.lammpstrj"),
        cut     = 15.0,
        pair    = ["model=BSpline,type=A:A,min=3.9,max=15.0,resolution=0.3,order=4",
                   "model=BSpline,type=A:B,min=3.9,max=15.0,resolution=0.3,order=4",
                   "model=BSpline,type=A:C,min=3.9,max=15.0,resolution=0.3,order=4",
                   "model=BSpline,type=B:B,min=3.9,max=15.0,resolution=0.3,order=4",
                   "model=BSpline,type=B:C,min=3.9,max=15.0,resolution=0.3,order=4",
                   "model=BSpline,type=C:C,min=3.9,max=15.0,resolution=0.3,order=4"
                  ],
        verbose = 0,
        save    = 'return'
    )

    benchmark = [0.351629, 0.991323, 1.969529, 3.101544, 3.419564, 3.550187, 3.641121, 3.800916, 4.143540, 4.428851, 4.725986, 5.158559, 5.641956, 6.197080, 6.909919, 7.641150, 8.623013, 9.450326,10.118991,10.727784,11.486072,12.172890,12.641396,13.293030,14.171136,15.210665,16.085361,16.772577,17.326560,17.991520,18.759520,19.491264,20.330532,21.120097,21.923910,22.729731,23.681647,18.352001,12.510831, 6.250291]
    
    coeffs = mean['Pair_A-A']
    print(','.join(["%9.6f" % (i) for i in coeffs]))
    print("")

    for i in range(len(benchmark)):
        diff = (coeffs[i] - benchmark[i]) / benchmark[i]
        print("X=%3d, Y0=%10.3e, Y=%10.3e, dY/Y0=%5.2f%%" %(i+1, benchmark[i], coeffs[i], diff*100))
        assert abs(diff)<0.01
    
    


"""
A-A
1.498569	4.193904	8.097605	12.416496	12.992205	12.905270	12.907778	13.128850	13.647911	14.555245	
15.793071	17.430897	19.388360	21.601030	24.352998	27.796083	31.429944	35.045197	37.704528	40.174972	
42.266003	44.561671	47.375890	50.348981	53.561326	56.886214	60.313042	63.622329	66.422979	69.646219	
72.715364	76.015658	79.651082	83.813377	87.560613	91.050251	95.378043	74.089301	50.676407	25.795335
"""