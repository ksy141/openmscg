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
      --model               Regularization settings (default: model=none)

'''

from mscg import *

def SolveNE(XtX, XtY, alpha = 0.0):
    if alpha > 1.0E-6:
        screen.info(["Solver => Ridge Regression"])
        
        scale = np.sqrt(np.square(XtX).sum(axis=1))
        XtX /= scale
        XtX = XtX + np.identity(XtX.shape[0]) * (alpha * alpha)
        
        #c = np.linalg.solve(XtX, XtY)
        c = np.matmul(np.linalg.pinv(XtX), XtY)
        c /= scale
    else:
        screen.info(["Solver => OLS"])
        #c = np.linalg.solve(XtX, XtY)
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
    group.add_argument("--alpha",  metavar='', type=float, default=0, help="alpha value for Ridge regression")
        
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
    
    c = SolveNE(X.copy(), y.copy(), args.alpha)
    screen.info(["Model coefficients:", c])
    
    if args.save == "return":
        return c
    
    else:
        offset = 0
        for name, m in models.items():
            m['params'] = c[offset:offset + m['params'].size]
            offset += m['params'].size
            
        Checkpoint(args.save, __file__).update({
            'models': models,
            'X': X,
            'y': y,
            'c': c
        }).dump()

    

if __name__ == '__main__':
    main()
