from mscg.cli import cgfm

cgfm.main(**{
    'top'    : "methanol_ucg.data",
    'traj'   : "cg.lammpstrj",
    'names'  : 'MeOH,Near,Far',
    'cut'    : 12.0,
    'pair'   : ['model=BSpline,type=Near:Near,min=3.0,max=10.0,resolution=0.1',
                'model=BSpline,type=Near:Far,min=3.0,max=10.0,resolution=0.1',
                'model=BSpline,type=Far:Far,min=3.0,max=10.0,resolution=0.1'],
    'ucg'    : 'replica=200,seed=1234',
    'ucg-wf' : 'RLE,target=MeOH,high=Near,low=Far,rth=4.5,wth=3.5',
    'verbose' : '1'
})
