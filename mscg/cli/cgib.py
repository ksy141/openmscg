
import numpy as np
from mscg import *


class Histogram:
    def __init__(self, n, name, args):
        self.ntype = n
        self.types = None
        self.name  = name
        
        self.min   = 0
        self.max   = 10
        self.bins  = 10
        
        self.id    = None
        self.x     = None
        self.n     = None
        self.U     = None
        
        segs = args.split(",")
        if len(segs) != self.ntype + 3:
            raise Exception('incorrect number of fields for option --' + self.name + ' ' + args)
        
        self.types = segs[:n]
        self.name = "-".join(self.types)
        
        for i in range(n, len(segs)):
            w = segs[i].split('=')
            if w[0] == "min":
                self.min = float(w[1])
            elif w[0] == "max":
                self.max = float(w[1])
            elif w[0] == "bins":
                self.bins = int(w[1])
            else:
                raise Exception('incorrect format of value for option --' + self.name + ' ' + segs[i])
        

def BuildHistAction(n, arg_name):
    class HistAction(argparse.Action):
        nbody = n
        name = arg_name
        
        def __call__(self, parser, namespace, values, option_string=None):
            
            getattr(namespace, self.dest).append(Histogram(HistAction.nbody, HistAction.name, values))
            return
        
        def help():
            msg = "define new " + HistAction.name + " analysis with format: "
            msg += ",".join(["type" + str(i+1) for i in range(HistAction.nbody)]) 
            msg += ",args; args and default values are: min=0,max=10,bins=10"
            return msg
            
    return HistAction



def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run MSCG Range-finder, RDF Calculations, or Inversed-boltzmann method. For detailed instructions please read ' + doc_root + 'commands/cgib.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--top",  metavar='file', action=TopAction, help=TopAction.help, required=True)
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--names",  metavar='', type=str, help="comma separated atom type names (needed when using LAMMPS data file for topology)")
    group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])
    group.add_argument("--cut",  metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    group.add_argument("--temp", metavar='', type=float, default=298.15, help="temperature (K) for IB")
    
    PairAction = BuildHistAction(2, "pair")
    group.add_argument("--pair", metavar='types,args', action=PairAction, help=PairAction.help(), default=[])
    
    BondAction = BuildHistAction(2, "bond")
    group.add_argument("--bond", metavar='types,args', action=BondAction, help=BondAction.help(), default=[])
    
    AngleAction = BuildHistAction(3, "angle")
    group.add_argument("--angle", metavar='types,args', action=AngleAction, help=AngleAction.help(), default=[])
    
    group.add_argument("--plot", metavar='U or N', type=str, default='U', help="plot the results of U (potential) or n (distribition)")
    
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
        
    # prepare plots

    if args.pair is not None:
        for pair in args.pair:
            screen.info("Add pair plot: " + pair.name)
            pair.id = args.top.get_pair_type(pair.types[0], pair.types[1])
            
    if args.bond is not None:
        for bond in args.bond:
            screen.info("Add bond plot: " + bond.name)
            bond.id = args.top.get_bond_type(bond.types[0], bond.types[1])

    if args.angle is not None:
        for angle in args.angle:
            screen.info("Add angle plot: " + angle.name)
            angles.id = args.top.get_angle_type(angle.types[0], angle.types[1], angle.types[2])

    # start processing trajectory
    
    TIMER.reset()
    last = TIMER.last
    
    for reader in args.traj:
        
        screen.info("Process trajectory: " + reader.file)
        trj = reader.traj
        
        if args.top.natoms != trj.natoms:
            screen.fatal("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (args.top.natoms, trj.natoms))
        
        cut2 = plist.cut * 2;
        if trj.box[0]<cut2 or trj.box[1]<cut2 or trj.box[2]<cut2:
            screen.fatal("Incorrect cut-off for the trajectory: cut-off (%f) must be larger than half of the box dimentions (%s)" % (plist.cut, str(trj.box)))
        
        plist.setup_bins(trj)
        start = TIMER.last
        
        while reader.next_frame():
            
            # process pair styles

            if len(args.pair)>0: 
                
                plist.build(trj)
                pstart = 0
                pn = 100000 # pairs per page
                
                while True:
                    z = plist.get_pairs(pstart, pn)
                    types = np.array(z[0])
                    vals = np.array(z[3])

                    for pair in args.pair:
                        dr = vals[types==pair.id]
                        hist, edges = np.histogram(dr, bins=pair.bins, range=(pair.min, pair.max))

                        if pair.n is None:
                            pair.n, pair.x = hist, edges[:-1] + np.diff(edges) * 0.5
                        else:
                            pair.n += hist

                    pstart += pn
                    if len(types)<pn:
                        break

            # process bonding styles

            if len(args.bond)>0 or len(args.angle)>0:
                blist.build(trj)

                def process_hist(one, types, vals):

                    vals = vals[types==one.id]
                    hist, edges = np.histogram(vals, bins=one.bins, range=(one.min, one.max))

                    if one.n is None:
                        one.n, one.x = hist, edges[:-1] + np.diff(edges) * 0.5
                    else:
                        one.n += hist

                for one in args.bond:
                    z = blist.get_bonds()
                    types = np.array(z[0])
                    dr = np.array(z[3])
                    process_hist(one, types, dr)

                for one in args.angle:
                    z = blist.get_angles()
                    types = np.array(z[0])
                    dr = np.array(z[4]) * (180.0 / np.pi)
                    process_hist(one, types, dr)
            
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

            # end of one frame
        
        # end of one trajectory
        
        if screen.verbose > 0:
            TIMER.click(None)
            elapsed = TIMER.last - start
            msg = " -> Processed %d frames. Elapsed: %0.0f secs." % (reader.nread, elapsed)
            print(('\r%s' % (msg)) + " " * 30)
        
    # end of processing trajectories
    
    # dump results
    
    def post_process(d, prefix):
        valid = d.n > 1.0E-40
        d.x = d.x[valid]
        d.n = d.n[valid]
        d.U = -0.0019872041 * args.temp *np.log(d.n)
        np.savetxt(prefix + '-' + d.name + '.dat', np.vstack([d.x, d.n, d.U]).T)
        return



    for pair in args.pair:
        pair.n = np.divide(pair.n, 4.0 * np.pi * np.square(pair.x))
        pair.n = np.divide(pair.n, pair.n[-1])
        post_process(pair, 'Pair')

    for bond in args.bond:
        bond.n = np.divide(bond.n, bond.n.max())
        post_process(bond, 'Bond')

    for angle in args.angle:
        angle.n = np.divide(angle.n, angle.n.max())
        post_process(angle, 'Angle')
    
    

    if args.plot != 'none':
        import matplotlib.pyplot as plt

        for pair in args.pair:
            plt.plot(pair.x, getattr(pair, args.plot), label='Pair ' + pair.name)
        
        for bond in args.bond:
            plt.plot(bond.x, getattr(bond, args.plot), label='Bond ' + bond.name)
        
        for angle in args.angle:
            plt.plot(angle.x, getattr(angle, args.plot), label='Angle ' + angle.name)
            
        plt.legend(loc='upper right')
        plt.show()
    
    screen.info("Processing is finished.")



if __name__ == '__main__':
    main()