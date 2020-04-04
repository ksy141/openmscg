
import numpy as np
from cg import *

def main(*args, **kwargs):
    
    # parse argument
    
    parser = CLIParser(description='Run MSCG Range-finder, RDF Calculations, or Inversed-boltzmann method.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--top",  metavar='file:format', type=str, help="topology file: in format of [filename,format:lammps|cg]", required=True)
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
    group.add_argument("--traj", metavar='', type=str, help="trajectory files: in format of [filename,skip:0,every:1,nread:0]", action='append')
    group.add_argument("--cut",    metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    group.add_argument("--temp",   metavar='', type=float, default=298.15, help="temperature (K) for IB")

    group.add_argument("--pair", metavar='', type=str, help="add a pair RDF, in format of [type1,type2,min,max,bins]", action='append')
    group.add_argument("--bond", metavar='', type=str, help="add a bond distribution, in format of [type1,type2,min,max,bins]", action='append')
    group.add_argument("--angle", metavar='', type=str, help="add an angle distribution, in format of [type1,type2,type3,min,max,bins]", action='append')
    group.add_argument("--plot", metavar='U|N', type=str, default='U', help="plot the results of U (potential) or n (distribition)")
    
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenCG CLI Command: " + __name__)
    
    # load topology
    
    top_args = [s.strip() for s in args.top.split(',')]
    
    if len(top_args)<2:
        screen.fatal(["wrong format of the topology argument.", "file format is not specified."])
    
    screen.info("Load topology: " + args.top)
    
    top = build_topology(top_args[1], top_args[0], {
        'names': [] if args.names is None else args.names.split(',')
    })

    screen.info("Generate bonds/angles/dihedrals ...")
    top.build_special(True, True, True)
    
    # prepare lists
    
    screen.info("Build pair and bonding list-based algorithm ...")
    plist = PairList(top)
    plist.init(cut = args.cut)
    blist = BondList(top)
        
    # prepare pair plots

    pairs = []

    if args.pair is not None:
        for pair in args.pair:
            screen.info("Add pair plot: " + pair)
            w = pair.split(',')
            pairs.append({
                'name': '-'.join(w[:2]),
                'type': top.get_pair_type(w[0], w[1]),
                'min' : float(w[2]),
                'max' : float(w[3]),
                'bins': int(w[4]),
                'n'   : None,
                'x'   : None,
            })

    # prepare bond plots

    bonds = []

    if args.bond is not None:
        for bond in args.bond:
            screen.info("Add bond plot: " + bond)
            w = bond.split(',')
            bonds.append({
                'name': '-'.join(w[:2]),
                'type': top.get_bond_type(w[0], w[1]),
                'min' : float(w[2]),
                'max' : float(w[3]),
                'bins': int(w[4]),
                'n'   : None,
                'x'   : None,
            })

    # prepare angle plots

    angles = []

    if args.angle is not None:
        for angle in args.angle:
            screen.info("Add angle plot: " + angle)
            w = angle.split(',')
            angles.append({
                'name': '-'.join(w[:3]),
                'type': top.get_angle_type(w[0], w[1], w[2]),
                'min' : float(w[3]),
                'max' : float(w[4]),
                'bins': int(w[5]),
                'n'   : None,
                'x'   : None,
            })

    # start processing trajectory
    
    TIMER.reset()
    last = TIMER.last
    
    for trj_arg in args.traj:
        
        screen.info("Process trajectory: " + trj_arg)
        
        w = trj_arg.split(',')
        filename, skip, every, frames = w[0], 0, 1, 0
        
        if len(w)>1:
            skip = int(w[1])
        
        if len(w)>2:
            every = int(w[2])
        
        if len(w)>3:
            frames = int(w[3])
                
        trj = Trajectory(filename)
        
        if top.natoms != trj.natoms:
            screen.fatal("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (top.natoms, trj.natoms))
        
        cut2 = plist.cut * 2;
        if trj.box[0]<cut2 or trj.box[1]<cut2 or trj.box[2]<cut2:
            screen.fatal("Incorrect cut-off for the trajectory: cut-off (%f) must be larger than half of the box dimentions (%s)" % (plist.cut, str(trj.box)))
        
        plist.setup_bins(trj)
        start = TIMER.last
        nread = 0
                
        for i in range(skip):
            trj.read_frame()
        
        while trj.read_frame():
            
            # process pair styles

            if len(pairs)>0: 
                
                plist.build(trj)
                pstart = 0
                pn = 100000 # pairs per page
                
                while True:
                    z = plist.get_pairs(pstart, pn)
                    types = np.array(z[0])
                    vals = np.array(z[3])

                    for pair in pairs:
                        dr = vals[types==pair['type']]
                        hist, edges = np.histogram(dr, bins=pair['bins'], range=(pair['min'], pair['max']))

                        if pair['n'] is None:
                            pair['n'], pair['x'] = hist, edges[:-1] + np.diff(edges) * 0.5
                        else:
                            pair['n'] += hist

                    pstart += pn
                    if len(types)<pn:
                        break

            # process bonding styles

            if len(bonds)>0 or len(angles)>0:
                blist.build(trj)

                def process_hist(one, types, vals):

                    vals = vals[types==one['type']]
                    hist, edges = np.histogram(vals, bins=one['bins'], range=(one['min'], one['max']))

                    if one['n'] is None:
                        one['n'], one['x'] = hist, edges[:-1] + np.diff(edges) * 0.5
                    else:
                        one['n'] += hist

                for one in bonds:
                    z = blist.get_bonds()
                    types = np.array(z[0])
                    dr = np.array(z[3])
                    process_hist(one, types, dr)

                for one in angles:
                    z = blist.get_angles()
                    types = np.array(z[0])
                    dr = np.array(z[4]) * (180.0 / np.pi)
                    process_hist(one, types, dr)

            # goto the next frame

            for i in range(every -1):
                trj.read_frame()

            nread += 1

            # timing

            if screen.verbose > 0:
                TIMER.click(None)

                if TIMER.last - last > 1.0:
                    last = TIMER.last
                    elapsed = TIMER.last - start

                    if frames>0:
                        remained = (TIMER.last - start) / nread * (frames - nread)
                        msg = " -> Processed %d of %d frames. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (nread, frames, elapsed, remained)
                    else:
                        msg = " -> Processed %d frames. Elapsed: %0.0f secs ..." % (nread, elapsed)

                    print('\r%s' % (msg), end="")

            # end of one frame

            if frames > 0 and nread >= frames:
                break
        
        # end of one trajectory
        
        if screen.verbose > 0:
            TIMER.click(None)
            elapsed = TIMER.last - start
            msg = " -> Processed %d frames. Elapsed: %0.0f secs." % (nread, elapsed)
            print(('\r%s' % (msg)) + " " * 30)
        
    # end of processing trajectories
    
    # dump results
    
    def post_process(d, prefix):
        valid = d['n'] > 1.0E-40
        d['x'] = d['x'][valid]
        d['n'] = d['n'][valid]
        d['U'] = -0.0019872041 * args.temp *np.log(d['n'])
        np.savetxt(prefix + '-' + d['name'] + '.dat', np.vstack([d['x'], d['n'], d['U']]).T)
        return



    for pair in pairs:
        pair['n'] = np.divide(pair['n'], 4.0 * np.pi * np.square(pair['x']))
        pair['n'] = np.divide(pair['n'], pair['n'][-1])
        post_process(pair, 'Pair')

    for bond in bonds:
        bond['n'] = np.divide(bond['n'], bond['n'].max())
        post_process(bond, 'Bond')

    for angle in angles:
        angle['n'] = np.divide(angle['n'], angle['n'].max())
        post_process(angle, 'Angle')



    if args.plot != 'none':
        import matplotlib.pyplot as plt

        for pair in pairs:
            plt.plot(pair['x'], pair[args.plot], label='Pair ' + pair['name'])
        
        for bond in bonds:
            plt.plot(bond['x'], bond[args.plot], label='Bond ' + bond['name'])
        
        for angle in angles:
            plt.plot(angle['x'], angle[args.plot], label='Bond ' + angle['name'])
            
        plt.legend(loc='upper right')
        plt.show()
    
    screen.info("Processing is finished.")



if __name__ == '__main__':
    main()