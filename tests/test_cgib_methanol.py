import pytest
from mscg.cli import cgib

def test_2s(datafile):
    
    dfs = cgib.main(
        top     = datafile("methanol_1728_2s.data"),
        traj    = datafile("methanol_1728_2s.trr,frames=500"),
        names   = 'CH3,OH',
        cut     = 10.0,
        pair    = ['CH3,CH3,min=2.8,max=10.0,bins=200', 'CH3,OH,min=2.8,max=10.0,bins=200', 'OH,OH,min=2.5,max=10.0,bins=200'],
        bond    = ['CH3,OH,min=1.35,max=1.65,bins=60'],
        verbose = 0,
        save    = 'return'
    )
        
    print("")
    
    benchmarks = [
        [-0.38171003, -0.38762647, -0.39988627, -0.40752347, -0.4089117,  -0.40874726, -0.4062567, -0.39776143, -0.39525766, -0.39181652],
        [-0.19976016, -0.17234927, -0.1486105,  -0.1224553,  -0.10083989, -0.0750504, -0.05117479, -0.03109205, -0.01499615,  0.00036485],
        [0.64377958, 0.62199957, 0.61026951, 0.57939914, 0.55971546, 0.52681372, 0.50739519, 0.4825158,  0.44840662, 0.42170519],
        [4.40205329e-04, 2.03246518e-02, 6.00544327e-02, 1.08639149e-01, 1.82976844e-01, 2.62163523e-01, 3.59980819e-01, 4.74700765e-01, 6.07324478e-01, 7.32104652e-01]
    ]
    
    for i in range(4):
        d = dfs[i]['data'][30:40,2]
        
        for j in range(10):
            assert abs(d[j]-benchmarks[i][j])<0.001