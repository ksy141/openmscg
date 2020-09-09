'''Calculate derivatives with respect to the model parameters from CG trajectories.

Description
-----------

*To be added*

Usage
-----

Syntax of running ``cgderiv`` command ::

    General arguments:
      -h, --help          show this help message and exit
      -v L, --verbose L   screen verbose level (default: 0)
      --save              file name for model output (default: model)

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
      --pair types,args   define new pair model with format:
                          style,type1,style,type2,kwargs (default: [])

'''

from mscg import *
import numpy as np

class CGDeriv:
    
    def __init__(self, *args, **kwargs):

        # parse argument

        desc = 'Calculate derivatives with respect to parameters in REM.' + doc_root + 'commands/cgderiv.html'

        parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)

        group = parser.add_argument_group('General arguments')
        group.add_argument("-h", "--help", action="help", help="show this help message and exit")
        group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")

        group.add_argument("--save",  metavar='', type=str, default="model", help="file name for model output")

        group = parser.add_argument_group('Required arguments')
        group.add_argument("--top",  metavar='file', action=TopAction, help="topology file", required=True)

        group = parser.add_argument_group('Optional arguments')
        group.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
        group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])

        group.add_argument("--cut", metavar='', type=float, default=10.0, help="cut-off for pair interactions")

        PairAction = BuildModelArgAction(2, "pair");
        group.add_argument("--pair", metavar='types,args', action=PairAction, help=PairAction.help(), default=[])
    
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

        screen.info("Generate bonds/angles/dihedrals ...")
        args.top.build_special(True, True, True)
        
        # prepare lists

        screen.info("Build pair and bonding list-based algorithm ...")
        plist = PairList(args.top)
        plist.init(cut = args.cut)
        blist = BondList(args.top)

        # build up tables
        model.empty()

        for pair in args.pair: pair.create(args.top, plist)
        for m in model.items: m.require('eval_deriv_u')
        
        # save references
        
        self.args = args
        self.plist = plist
        self.blist = blist
        
        
    
    def process(self):
        
        # init dudl
        dudl_frames = {}
        
        for m in model.items:
            dudl_frames[m.name] = []
        
        # start processing trajectory
        
        TIMER.reset()
        last = TIMER.last
        
        for reader in TrajBatch(self.args.traj, natoms = self.args.top.natoms, cut = self.plist.cut):

            if reader.nread == 1:
                self.plist.setup_bins(reader.traj)

            TIMER.click('io')
            TIMER.click('pair', self.plist.build(reader.traj))
            TIMER.click('bond', self.blist.build(reader.traj))
            TIMER.click('eval', model.call('eval_deriv_u'))
            
            for m in model.items:
                dudl_frames[m.name].append(m.dudl.copy())

        # end of processing trajectories

        screen.info([""] + TIMER.report(False) + [""])

        # end
        
        dudl_mean, dudl_var = {}, {}
        
        for m in model.items:
            if len(dudl_frames[m.name]) == 0:
                continue
            
            dudls = np.stack(dudl_frames[m.name], axis=0);
            dudl_mean[m.name] = dudls.mean(axis=0)
            dudl_var[m.name] = dudls.var(axis=0)
            
        # collect results

        if self.args.save == 'return':
            return dudl_mean, dudl_var
        else:
            Checkpoint(self.args.save).update({
                'models'    : model.serialize(),
                'dudl_mean' : dudl_mean,
                'dudl_var'  : dudl_var
            }).dump()
            
            
        
def main(*args, **kwargs):
    return CGDeriv(*args, **kwargs).process()

if __name__ == '__main__':
    main()
