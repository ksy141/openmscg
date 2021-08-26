'''Normal Equations Solver for Force-Matching Method

Description
-----------

The ``cgnes`` provides the module to solve the normal equations problem from 
force-matching approach. It can be also used in processing the stored normal
equations from previous FM tasks, under different conditions, such as the
Alpha value for Ridge regularization. Another important usage of this tool
is to read multile saved equations for the same modeling task but constructed
from different parts of trajectories, make an elemetary summation of them, and
solve the NE problem of it. This is the last step for the `distributed workflow
<../performance.html#distributed-workflow>`__ for FM/UCG modeling. 

Usage
-----

Syntax of running ``cgnes`` command ::

    usage: cgnes.py [-h] [-v L] [--save] [--equation files [files ...]] [--alpha]
                    [--bayesian]

    General arguments:
      -h, --help            show this help message and exit
      -v L, --verbose L     screen verbose level (default: 0)

    Optional arguments:
      --save                file name for matrix output (default: result)
      --equation files [files ...]
                            CGFM result files (default: None)
      --alpha               alpha value for Ridge regression (default: 0)
      --bayesian            maximum steps for Bayesian regularization (default: 0)

'''

from mscg import *

def SolveNE(nmmat, alpha = 0.0, bayesian = 0):
    XtX, XtY = nmmat['XtX'].copy(), nmmat['XtY'].copy()
    
    def GetR(XtX, XtY, c):
        return np.dot(c, np.matmul(XtX, c)) - 2.0 * np.dot(XtY, c) + nmmat['y_sumsq']
        
    def Normalize(X):
        scale = np.sqrt(np.square(X).sum(axis=1))
        return X / scale, scale
        
    if alpha > 1.0E-6:
        screen.info(["Solver => Ridge Regression"])
        
        XtX_, scale = Normalize(XtX)
        XtX_ = XtX_ + np.identity(XtX_.shape[0]) * (alpha * alpha)
        
        c = np.matmul(np.linalg.pinv(XtX_), XtY) / scale
    else:
        screen.info(["Solver => OLS"])
        c = np.matmul(np.linalg.pinv(XtX), XtY)
        
    if bayesian > 0:
        screen.info(["Applying Bayesian iterative regularization ..."])
        N = nmmat['row_per_frame'] * nmmat['nframe']
        
        _a = XtX.shape[0] / np.dot(c, c)
        _b = N / GetR(XtX, XtY, c)
        a = np.ones(XtX.shape[0]) * (_a)
        screen.info(["Bayesian iteration step %d: beta = %lf" % (0, _b)])
        
        
        for it in range(bayesian):
            
            XtX_reg = XtX + np.identity(XtX.shape[0]) * (a / _b)
            XtX_norm, scale = Normalize(XtX_reg)
            c = np.matmul(np.linalg.pinv(XtX_norm), XtY) / scale 
            
            invXtX_reg = np.linalg.pinv(XtX_reg)
            a = 1.0 / (c * c + np.diagonal(invXtX_reg)/_b)
            
            tr = np.trace(np.matmul(invXtX_reg, XtX))
            _b = (N - tr) / GetR(XtX, XtY, c)
            screen.info(["Bayesian iteration step %d: beta = %lf" % (it+1, _b)])
            
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
    group.add_argument("--bayesian",  metavar='', type=int, default=0, help="maximum steps for Bayesian regularization")
    
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
    matrix = chk.data['matrix']
    
    for eq in args.equation[1:]:
        data = Checkpoint.load(eq).data
        matrix['XtX'] += data['matrix']['XtX']
        matrix['XtY'] += data['matrix']['XtY']
        matrix['y_sumsq'] += data['matrix']['y_sumsq']
        matrix['nframe'] += data['matrix']['nframe']
    
    c = SolveNE(matrix, args.alpha, args.bayesian)
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
            'matrix': matrix,
            'c': c
        }).dump()

    

if __name__ == '__main__':
    main()
