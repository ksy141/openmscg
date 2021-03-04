'''Normal Equations Solver for Force-Matching Method

Description
-----------

The ``cgnes`` provides the module to solve the normal equations problem from 
force-matching approach. It can be also used in processing the stored normal
equations from previous FM tasks, under different conditions, such as the
Lambda value for LASSO regularization. Another important usage of this tool
is to read multile saved equations for the same modeling task but constructed
from different parts of trajectories, make an elemetary summation of them, and
solve the NE problem of it. This is the last step for the `distributed workflow
<../performance.html#distributed-workflow>`__ for FM/UCG modeling. 

Usage
-----

Syntax of running ``cgnes`` command ::

    usage: cgnes.py [-h] [-v L] [--save] [--equation files [files ...]] [--lasso]

    General arguments:
      -h, --help            show this help message and exit
      -v L, --verbose L     screen verbose level (default: 0)

    Optional arguments:
      --save                file name for matrix output (default: result)
      --equation files [files ...]
                            CGFM result files (default: None)
      --lasso               lambda value for Lasso regularizer (default: 0.0)

'''

from mscg import *

def SolveNE(XtX, XtY, Lasso=None):
    
    if Lasso > 1.0e-32:
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
    
    return c


def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run normal equations solver for FM/UCG method. For detailed instructions please read ' + doc_root + 'commands/cgnes.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--save",  metavar='', type=str, default="result", help="file name for matrix output")
    group.add_argument("--equation", metavar='files', nargs='+', type=argparse.FileType('rb'), help="CGFM result files")
    group.add_argument("--lasso",  metavar='', type=float, default=0.0, help="lambda value for Lasso regularizer")
        
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)
    
    # load NEs
    
    if args.equation is None:
        raise Exception('No FM result file with normal equations is provided.')
    
    chk = Checkpoint.load(args.equation[0])
    models = chk.data['models']
    X = chk.data['X']
    y = chk.data['y']
    
    for eq in args.equation[1:]:
        data = Checkpoint.load(eq).data
        X += data['X']
        y += data['y']
    
    c = SolveNE(X, y, args.lasso)
    screen.info(["Model coefficients:", c])
    
    if args.save == "return":
        return c
    
    else:
        Checkpoint(args.save, __file__).update({
            'models': models,
            'X': X,
            'y': y,
            'c': c
        }).dump()

    

if __name__ == '__main__':
    main()
