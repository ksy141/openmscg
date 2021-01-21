'''Force-matching to construct force tables with B-Spline or other linear models

Description
-----------

The ``cgfm`` command is used to do the force-matching method with a given CG topology and a set of CG trajectories. It is the main computing engine of the OpenMSCG package. Briefly, it reads in the coordinates and forces from the trajectories, calculating the loadings for a group of user-defined force tables, storing them in a covariance matrix, and finaly solving the linear regression to get the coefficients for each force types. Based on the variational theory, the computed coefficients defines the force function that minimize the force residues between real CG values and models. 

In OpenMSCG, a force table is defined in B-spline function defined by *n* **uniformly** spaced knots. Users can specify the starting value (``min``), ending value (``max``), the interval between adjacent knonts (``resolution``) and the ``order`` of the spline. Finally, the number of knots, which is also the number of computed coefficients, is

``n = (max - min) / resolution + order - 2``

Usage
-----

Syntax of running ``cgfm`` command ::

    usage: cgfm [-h] [-v L] --top file [--names] [--traj file[,args]] [--cut]
                [--save] [--lasso] [--pair [key=value]] [--bond [key=value]]
                [--angle [key=value]] [--dihedral [key=value]] [--ucg [key=value]]
                [--ucg-wf [key=value]]

    Run MSCG force-matching method. For detailed instructions please read
    https://software.rcc.uchicago.edu/mscg/docs/commands/cgfm.html

    General arguments:
      -h, --help            show this help message and exit
      -v L, --verbose L     screen verbose level (default: 0)

    Required arguments:
      --top file            topology file (default: None)

    Optional arguments:
      --names               comma separated atom names (needed when using LAMMPS
                            data file for topology) (default: None)
      --traj file[,args]    reader for a trajectory file, multiple fields
                            separated by commas, the first field is the file name,
                            while others define the skip, every and frames
                            (default args: file,skip=0,every=1,frames=0) (default:
                            [])
      --cut                 cut-off for pair interactions (default: 10.0)
      --save                file name for matrix output (default: result)
      --lasso               lambda value for Lasso regularizer (default: 0.0)
      --pair [key=value]    add a model declaration for pair-style interactions.
                            (default: [])
      --bond [key=value]    add a model declaration for bond-style interactions.
                            (default: [])
      --angle [key=value]   add a model declaration for angle-style interactions.
                            (default: [])
      --dihedral [key=value]
                            add a model declaration for dihedral-style
                            interactions. (default: [])
      --ucg [key=value]     settings for UCG modeling (default: None)
      --ucg-wf [key=value]  define new state-function for UCG (default: [])

'''

from mscg import *

def main(*args, **kwargs):
    
    models.empty()
    
    # parse argument
    
    desc = 'Run MSCG force-matching method. For detailed instructions please read ' + doc_root + 'commands/cgfm.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--top",  metavar='file', action=TopAction, help="topology file", required=True)
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
    group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])
    
    group.add_argument("--cut", metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    group.add_argument("--save",  metavar='', type=str, default="result", help="file name for matrix output")
    
    group.add_argument("--lasso",  metavar='', type=float, default=0.0, help="lambda value for Lasso regularizer")
    
    group.add_argument("--pair",  metavar='[key=value]', action=ModelArgAction, help=ModelArgAction.help('pair'), default=[])
    group.add_argument("--bond",  metavar='[key=value]', action=ModelArgAction, help=ModelArgAction.help('bond'), default=[])
    group.add_argument("--angle", metavar='[key=value]', action=ModelArgAction, help=ModelArgAction.help('angle'), default=[])
    group.add_argument("--dihedral", metavar='[key=value]', action=ModelArgAction, help=ModelArgAction.help('dihedral'), default=[])
        
    group.add_argument("--ucg", metavar='[key=value]', action=UCGArgAction, help=UCGArgAction.help(), default=None)
    group.add_argument("--ucg-wf", metavar='[key=value]', action=WFArgAction, help=WFArgAction.help(), default=[])
    
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenCG CLI Command: " + __name__)
    
    # load topology
    
    screen.info("Check topology ... ")
    
    if args.names is not None:
        args.top.reset_names(args.names.split(','))
    
    # prepare lists
    
    screen.info("Build pair and bonding list-based algorithm ...")
    plist = PairList(cut = args.cut, binsize = args.cut * 0.5)
    plist.init(args.top.types_atom, args.top.linking_map(True, True, True))
    blist = BondList(
        args.top.types_bond, args.top.bond_atoms, 
        args.top.types_angle, args.top.angle_atoms, 
        args.top.types_dihedral, args.top.dihedral_atoms)
    UCG.init(plist, blist)
    
    # build up tables
        
    [pair.setup(args.top, plist) for pair in args.pair]
    [bond.setup(args.top, blist) for bond in args.bond]
    [angle.setup(args.top, blist) for angle in args.angle]
    [dihedral.setup(args.top, blist) for dihedral in args.dihedral]
    
    for model in models.items:
        screen.info(" ".join([str(i) for i in 
            ["Model:", model.style, model.name, "T-" + str(model.tid)]]))
    
    # build up coefficients matrix
    
    screen.info("Build coefficients matrix ...")
    
    n = sum([model.nparam for model in models.items])
    matrix_cov = np.zeros(shape=(n, n))
    vector_cov = np.zeros(n)
    
    # start processing trajectory
    
    TIMER.reset()
    last = TIMER.last
            
    for reader in TrajBatch(args.traj, natoms = args.top.n_atom, cut = plist.cut):
    
        if reader.nread == 1:
            plist.setup_bins(reader.traj.box)
        
        TIMER.click('io')
        TIMER.click('pair', plist.build(reader.traj.x))
        TIMER.click('bond', blist.build(reader.traj.box, reader.traj.x))
        
        for types in UCGSpawner(args.top, args.traj):
            if types is not None:
                plist.update_types(np.array(types, dtype=np.int32))
                TIMER.click('ucg')
            
            TIMER.click('model', models.compute_fm())
            
            matrix_coeff = np.hstack([model.dF for model in models.items])
            matrix_cov += np.matmul(matrix_coeff.T, matrix_coeff)
            vector_cov += np.matmul(matrix_coeff.T, reader.traj.f.astype(np.float64).flatten())
            TIMER.click('matrix')
        
    # end of processing trajectories
    
    XtX = matrix_cov
    XtY = vector_cov
    
    if args.lasso > 1.0e-32:
        screen.info(["Solver => LASSO_LARS"])
        
        try:
            from sklearn import linear_model
            clf = linear_model.LassoLars(alpha=args.lasso, fit_intercept=False, max_iter=50000)
        except:
            screen.fatal("Package [sklearn] is required when using the LASSO estimator. (See https://scikit-learn.org/stable/install.html)")
            
        clf.fit(XtX, XtY)
        c = clf.coef_
    else:
        screen.info(["Solver => OLS"])
        c = np.matmul(np.linalg.pinv(XtX), XtY)
    
    screen.info(["Model coefficients:", c])
    
    TIMER.click('solver')
        
    if args.save != "return":
        offset = 0
        
        for m in models.items:
            m.params[:] = c[offset:offset + m.nparam]
            offset += m.nparam
        
        Checkpoint(args.save, __file__).update({'models': models.serialize(), 'X': matrix_cov, 'c': c}).dump()
    
    screen.info([""] + TIMER.report(False) + [""])
    
    # end
    
    if args.save == "return":
        return c
    

if __name__ == '__main__':
    main()
