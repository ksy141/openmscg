from mscg.cli import cgfm

cgfm.main(
    top  = "cg.top",
    traj = "cg.lammpstrj",
    cut  = 10.0,
    pair = ['model=BSpline,type=MeOH:MeOH,min=2.8,max=7.0,resolution=0.1,order=3'],
    verbose = '1'
)
