'''Iterative relative-entroy method (REM) for CG modeling.

Description
-----------

The ``cgrem`` command the main engine that drives the REM approach. It applies the iterative approch to minimize the differences of ensemble averages of between derivatives from the reference model and trial models. Therefore, this command calls the ``cgderiv`` command and an MD engine (currently only LAMMPS is supported) in loops.

Usage
-----

Syntax of running ``cgrem`` command ::

    General arguments:
      -h, --help            show this help message and exit
      -v L, --verbose L     screen verbose level (default: 0)

    Required arguments:
      --ref file            checkpoint file for reference model (default: None)
      --cgderiv-arg file    file for cgderiv arguments (default: None)
      --md file             file containing MD command (default: None)

    Optional arguments:
      --chi                 Chi-value (default: 1.0)
      --table               prefix of table names (default: )
      --maxiter             maximum iterations (default: 20)
      --restart file        restart file (default: restart)
      --models file         initial model parameters (default: model.txt)
      --optimizer name,args Define optimizer (default: [])

'''

from mscg import *
from mscg.cli.cgderiv import CGDeriv
from mscg.table import Table

import pickle, os, sys

class OptimizerBuiltin:
    
    def __init__(self, **kwargs):
        
        self.chi = float(kwargs.get('chi', 0.05))
        self.t = float(kwargs.get('t', 298.15))
        self.beta = 1.0 / (0.001985875 * self.t)        
        
    def run(self, params, dudl_ref, dudl_mean, dudl_var):
        if dudl_mean is None:
            return params

        for name, dudl_aa in dudl_ref.items():
            dudl_cg = dudl_mean[name].copy()
            var = dudl_var[name].copy()
            
            step = self.chi * (dudl_cg - dudl_aa) / (self.beta * var)
            param_prev = params[name].copy()
            params[name] = param_prev + step
            
            screen.info("")
            screen.info("Beta: " + str(self.beta))
            screen.info("chi: " + str(self.chi))
            
            screen.info("\n=> %15s %15s %15s" % ("<dU/dL>_aa", "<dU/dL>_cg", "Var(dU/dL)_cg"))
            for i in range(params[name].shape[0]):
                screen.info("=> %15.5e %15.5e %15.5e" % (dudl_aa[i], dudl_cg[i], var[i]))
            
            screen.info("\n=> %15s %15s %15s" % ("param_pre", "param_chg", "param_cur"))
            for i in range(params[name].shape[0]):
                screen.info("=> %15.5e %15.5e %15.5e" % (param_prev[i], step[i], params[name][i]))
            
        return params


class OptimizerAction(argparse.Action):
    
    help = "Name and parameters for the optimizer"
    
    def __call__(self, parser, namespace, values, option_string=None):
        
        if type(values) != str:
            raise ValueError("incorrect format of value for option --optimizer")
        
        segs = values.split(",")
        kwargs = {}
        
        for i in range(1, len(segs)):
            w = [seg.strip() for seg in segs[i].split("=")]
            
            if len(w)==1:
                raise ValueError("incorrect format of value for option --optimizer: " + segs[i])
            else:
                kwargs[w[0]] = w[1]
        
        if segs[0]=='builtin':
            model_class = OptimizerBuiltin
        else:
            sys.path.append(os. getcwd())
            model_module = importlib.import_module(segs[0])
            model_class = getattr(model_module, "Optimizer")
            
        setattr(namespace, self.dest, model_class(**kwargs))


def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Iterative REM. ' + doc_root + 'commands/cgrem.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--ref",  metavar='file', help="checkpoint file for reference model", required=True)
    group.add_argument("--cgderiv-arg", metavar='file', help="file for cgderiv arguments", required=True)
    group.add_argument("--md",  metavar='file', type=str, help="file containing MD command", required=True)
        
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--table", metavar='', type=str, default='', help="prefix of table names")
    group.add_argument("--maxiter", metavar='', type=int, default=20, help="maximum iterations")
    
    group.add_argument("--restart", metavar='file', default="restart", help="restart file")
    group.add_argument("--models", metavar='file', default="model.txt", help="initial model parameters")
    group.add_argument("--optimizer", metavar='name,[key=value]', action=OptimizerAction, help=OptimizerAction.help, default=[])
    
    # parse args
    
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)
    
    # build cgderiv object
    
    screen.info("Build up calculator for derivatives ...")
    deriv = CGDeriv('@' + args.cgderiv_arg)
    deriv.args.verbose = args.verbose
    deriv.args.save = 'return'
    screen.verbose = args.verbose
    
    # read reference model
    
    screen.info("Read reference model ...")
    ref = pickle.load(open(args.ref, 'rb'))
    
    for name, m in ref['models'].items():
        screen.info(" -> Model: " + name)
        
        if models[name] is None:
            screen.fatal("Reference model [%s] is not in targets." % (name))
            
        if m['nparam'] != models[name].nparam:
            screen.fatal("Incorrect number of parameters for reference model %s. (%d != %d)" % (name, models[name].nparam, m['nparam']))
    
    dudl_ref = ref['dudl_mean']
    
    # read tables
    
    targets = {}
    
    with open(args.models, "r") as f:
        rows = [_.strip().split() for _ in f.read().strip().split("\n")]

        for row in rows:
            model_name = row[0]
            model_params = [float(_) for _ in row[5:]]
            m = models[model_name]

            if m is None:
                screen.fatal("Model %s is not defined." % (model_name))

            if m.nparam != len(model_params):
                screen.fatal("Incorrect number of parameters for model %s. (%d != %d)" % (model_name, m.nparam, len(model_params)))
            
            if ('U' not in row[4]) and ('L' not in row[4]) and ('H' not in row[4]):
                screen.fatal("Incorrect padding option model %s. (%s)" % (model_name, row[4]))
            
            targets[model_name] = {
                'min': float(row[1]),
                'max': float(row[2]),
                'inc': float(row[3]),
                'pad': row[4],
                'init_params': model_params
            }
    
    # read md command
    
    with open(args.md, "r") as f:
        md_cmd = f.readline().strip()
        screen.info("MD command: " + md_cmd)
    
    # restart or init
    
    if os.path.isfile(args.restart + ".p"):
        screen.info("Restart iterations from checkpoint ...")
        iters     = pickle.load(open(args.restart + ".p", 'rb'))['iterations']
        params    = iters[-1]['params']
        dudl_mean = iters[-1]['dudl_mean']
        dudl_var  = iters[-1]['dudl_var']
    else:
        screen.info("Initialize models for iterations ...")
        iters    = []
        params   = {}
        dudl_mean, dudl_var = None, None
        
        for m in models.items:
            if m.name not in targets:
                screen.fatal("Parameters are not initialized for model %s" % (model_name))
        
            params[m.name] = targets[m.name]['init_params']
            
    # iterations
    
    for it in range(args.maxiter):
        screen.info("* Interation %d" % (len(iters)))
        
        # generate tables
        
        Checkpoint(args.restart + ".bak", __file__).update({'dudl_ref': dudl_ref, 'iterations': iters}).dump()
        params = args.optimizer.run(params.copy(), dudl_ref, dudl_mean, dudl_var)
                        
        for m in models.items:
            m.params = np.array(params[m.name])
            screen.info("Generate table [%s] to [%s]" % (m.name, args.table + m.name + '.table'))
            
            tbl = Table(m, force=False, prefix=args.table)
            tbl.compute(targets[m.name]['min'], targets[m.name]['max'], targets[m.name]['inc'])
            
            pad = targets[m.name]['pad']
            if "L2" in pad: tbl.padding_low2(targets[m.name]['min'])
            elif "L" in pad: tbl.padding_low(targets[m.name]['min'])
            if "H" in pad: tbl.padding_high(targets[m.name]['max'])
            
            tbl.dump_lammps()

        # run MD
        
        if os.system(md_cmd) != 0:
            screen.fatal("MD command terminated with failures.")
        
        # calculate derivative
        
        dudl_mean, dudl_var = deriv.process()
        
        # update checkpoint 
        
        iters.append({
            'params': params,
            'dudl_mean': dudl_mean,
            'dudl_var': dudl_var,
        })
        
        Checkpoint(args.restart, __file__).update({ 'dudl_ref': dudl_ref, 'iterations': iters}).dump()

    

if __name__ == '__main__':
    main()
