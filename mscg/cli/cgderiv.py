'''Calculate derivatives with respect to the model parameters from CG trajectories.

Description
-----------

The ``cgderiv`` command is used to calculate the derivates of potential energies on model parameters (*Lambdas*) from give trajectories, as well as corresponding variances. This is the key step for the relative-entroy method (REM), because REM is aimed to minimize the differences the means of such derivates between the reference trajectory (usually from the all-atom simulations) and trial trajectory from the models.

Usage
-----

Syntax of running ``cgderiv`` command ::

    usage: cgderiv.py [-h] [-v L] [--save] --top file [--names]
                      [--traj file[,args]] [--cut] [--pair [key=value]]

    General arguments:
      -h, --help          show this help message and exit
      -v L, --verbose L   screen verbose level (default: 0)

    Required arguments:
      --top file          topology file (default: None)

    Optional arguments:
      --names             comma separated atom names (needed when using LAMMPS
                          data file for topology) (default: None)
      --traj file[,args]  reader for a trajectory file, multiple fields separated
                          by commas, the first field is the file name, while
                          others define the skip, every and frames (default args:
                          file,skip=0,every=1,frames=0) (default: [])
      --cut               cut-off for pair interactions (default: 10.0)
      --pair [key=value]  add a model declaration for pair-style interactions.
                          (default: [])
      --save              file name for model output (default: model)

**Note**: modeling of bonded interactions hasn't been implemented yet.

'''

from mscg import *
import numpy as np

class CGDeriv:
    
    def __init__(self, *args, **kwargs):
        
        models.empty()
        
        # parse argument

        desc = 'Calculate derivatives with respect to parameters in REM.' + doc_root + 'commands/cgderiv.html'

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

        group.add_argument("--pair",  metavar='[key=value]', action=ModelArgAction, help=ModelArgAction.help('pair'), default=[])
        
        group.add_argument("--save",  metavar='', type=str, default="model", help="file name for model output")
    
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
        
        # setup models
        
        [pair.setup(args.top, plist) for pair in args.pair]
        
        # save references
        
        self.args = args
        self.plist = plist
        #self.blist = blist
            
    def process(self):
        
        # init dudl
        dudl_frames = {}
        
        for m in models.items:
            dudl_frames[m.name] = []
        
        # start processing trajectory
        
        TIMER.reset()
        last = TIMER.last
        
        for reader in TrajBatch(self.args.traj, natoms = self.args.top.n_atom, cut = self.plist.cut):
            
            if reader.nread == 1:
                self.plist.setup_bins(reader.traj.box)
            
            TIMER.click('io')
            TIMER.click('pair', self.plist.build(reader.traj.x))
            #TIMER.click('bond', blist.build(reader.traj.box, reader.traj.x))
            TIMER.click('model', models.compute_rem())
            
            for m in models.items:
                dudl_frames[m.name].append(m.dU.copy())

        # end of processing trajectories
        
        dudl_mean, dudl_var = {}, {}
        
        for m in models.items:
            dudls = np.stack(dudl_frames[m.name], axis=0);
            dudl_mean[m.name] = dudls.mean(axis=0)
            dudl_var[m.name] = dudls.var(axis=0)
        
        # collect results
        
        if self.args.save == 'return':
            return dudl_mean, dudl_var
        else:
            Checkpoint(self.args.save, __file__).update({
                'models'    : models.serialize(),
                'dudl_mean' : dudl_mean,
                'dudl_var'  : dudl_var
            }).dump()
        
        # end
        screen.info([""] + TIMER.report(False) + [""])
            
            
        
def main(*args, **kwargs):
    return CGDeriv(*args, **kwargs).process()

if __name__ == '__main__':
    main()
