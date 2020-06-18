
from mscg import *


class DeltaCreator:
    def __init__(self, n, style, args):
        self.ntype = n
        self.style = style
        
        segs = args.split(",")
        
        if len(segs) < self.ntype + 1:
            raise Exception('incorrect number of fields for option --' + self.style + ' ' + args)
        
        self.func = segs[0]
        self.types = segs[1:n+1]
        self.name = "_".join(self.types)
        
        self.kwargs = {}
        
        for i in range(n+1, len(segs)):
            w = segs[i].split('=')
            self.kwargs[w[0]] = w[1]
    
    def create(self, top, vlist):
        screen.info("Add %s coefficients delta: %s" % (self.style, self.name))
        
        get_type = getattr(top, "get_" + self.style + "_type")
        args = (vlist, self.name, get_type(*(self.types)))
        
        delta = globals()["Delta" + self.style.capitalize() + self.func]
        return delta(*args, **(self.kwargs)).setup_cache()



def BuildDeltaAction(n, arg_name):
    class DeltaAction(argparse.Action):
        nbody = n
        name = arg_name
        
        def __call__(self, parser, namespace, values, option_string=None):
            
            getattr(namespace, self.dest).append(DeltaCreator(DeltaAction.nbody, DeltaAction.name, values))
            return
        
        def help():
            msg = "define new " + DeltaAction.name + " model for delta_U/delta_eta"
            return msg
            
    return DeltaAction



def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run REM method. For detailed instructions please read ' + doc_root + 'commands/cgdev.html'
    
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
    
    PairAction = BuildDeltaAction(2, "pair");
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
    #blist = BondList(args.top)
    
    # build up tables
    
    pairs = [pair.create(args.top, plist) for pair in args.pair]

    # start processing trajectory
    
    TIMER.reset()
    last = TIMER.last
            
    for reader in args.traj:
        
        screen.info("Process trajectory: " + reader.file)
        trj = reader.traj
        
        if args.top.natoms != trj.natoms:
            screen.fatal("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (top.natoms, trj.natoms))
        
        cut2 = plist.cut * 2;
        if trj.box[0]<cut2 or trj.box[1]<cut2 or trj.box[2]<cut2:
            screen.fatal("Incorrect cut-off for the trajectory: cut-off (%f) must be larger than half of the box dimentions (%s)" % (plist.cut, str(trj.box)))
        
        plist.setup_bins(trj)
        start = TIMER.last

        while reader.next_frame():
            TIMER.click('io')
            TIMER.click('pair', plist.build(trj))
            #TIMER.click('bond', blist.build(trj))
            
            
            # timing

            if screen.verbose > 0:
                TIMER.click(None)

                if TIMER.last - last > 1.0:
                    last = TIMER.last
                    elapsed = TIMER.last - start

                    if reader.frames>0:
                        remained = (TIMER.last - start) / reader.nread * (reader.frames - reader.nread)
                        msg = " -> Processed %d of %d frames. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (reader.nread, reader.frames, elapsed, remained)
                    else:
                        msg = " -> Processed %d frames. Elapsed: %0.0f secs ..." % (reader.nread, elapsed)

                    print('\r%s' % (msg), end="")
        
        # end of one trajectory
        
        if screen.verbose > 0:
            TIMER.click(None)
            elapsed = TIMER.last - start
            msg = " -> Processed %d frames. Elapsed: %0.0f secs." % (reader.nread, elapsed)
            print(('\r%s' % (msg)) + " " * 30)
        
    # end of processing trajectories
        
    if args.save != "return":
        pass
    
    screen.info([""] + TIMER.report(False) + [""])
    
    # end
    
    if args.save == "return":
        pass
    

if __name__ == '__main__':
    main()
