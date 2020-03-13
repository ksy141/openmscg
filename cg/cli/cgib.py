
import numpy as np
from cg import *

import argparse

class CLIParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        if arg_line.strip()[:1] == '#':
            return []
        return arg_line.strip().split()

def main():
    
    parser = CLIParser(description='Run MSCG Range-finder, RDF Calculations, or Inversed-boltzmann method.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@')

    parser.add_argument("topology",   help="topology file",   type=str)
    parser.add_argument("trajectory", help="trajectory file", type=str)

    parser.add_argument("--skip",   metavar='N', type=int, default=0, help="skip first N frames")
    parser.add_argument("--every",  metavar='N', type=int, default=1, help="read the first frame in every N frames")
    parser.add_argument("--frames", metavar='N', type=int, default=0, help="maximum frames to be processed, 0 for all")

    parser.add_argument("--top",    metavar='', type=str, default="cgtop", help="type of topology file")
    parser.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
    parser.add_argument("--cut",    metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    parser.add_argument("--temp",   metavar='', type=float, default=298.15, help="temperature (K) for IB")

    parser.add_argument("--pair", metavar='', type=str, help="add a pair RDF, in format of [type1,type2,min,max,bins]", action='append')
    parser.add_argument("--bond", metavar='', type=str, help="add a bond distribution, in format of [type1,type2,min,max,bins]", action='append')
    parser.add_argument("--angle", metavar='', type=str, help="add an angle distribution, in format of [type1,type2,type3,min,max,bins]", action='append')
    parser.add_argument("--plot", metavar='', type=str, default='U', help="plot the results of U (potential) or n (distribition)")
    parser.add_argument("-v", "--verbose", help="Verbose mode", action="store_true")
    
    args = parser.parse_args()

    top = build_topology(args.top, args.topology, {
        'names': [] if args.names is None else args.names.split(',')
    })

    top.build_special(True, True, True)
    trj = Trajectory(args.trajectory)
    
    if top.natoms != trj.natoms:
        raise Exception("Inconsistent number of atoms between topology (%d) and trajectory (%d)." % (top.natoms, trj.natoms))

    plist = PairList(top, trj)
    plist.init(cut=args.cut)
    plist.setup_bins()
    
    blist = BondList(top, trj)
    nread = 0
    
    # prepare pair plots

    pairs = []

    if args.pair is not None:
        for pair in args.pair:
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

    for i in range(args.skip):
        trj.read_frame()
    
    if args.verbose:
        TIMER.reset()
        last = TIMER.last

    while trj.read_frame():
        
        # process pair styles

        if len(pairs)>0: 
            plist.build()
            pstart = 0
            pn = 100000

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
            blist.build()

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

        for i in range(args.every -1):
            trj.read_frame()

        nread += 1
        
        # timing
        
        if args.verbose:
            TIMER.click(None)

            if TIMER.last - last > 1.0:
                last = TIMER.last
                elapsed = TIMER.last - TIMER.start

                if args.frames>0:
                    remained = (TIMER.last - TIMER.start) / nread * (args.frames - nread)
                    msg = "Processed %d of %d frames. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (nread, args.frames, elapsed, remained)
                else:
                    msg = "Processed %d frames. Elapsed: %0.0f secs ..." % (nread, elapsed)

                print('\r%s' % (msg), end="")

        # end
        
        if args.frames>0 and nread >= args.frames:
            break
    
    print()


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



    if args.plot != '':
        import matplotlib.pyplot as plt

        for pair in pairs:
            plt.plot(pair['x'], pair[args.plot], label='Pair ' + pair['name'])
        
        for bond in bonds:
            plt.plot(bond['x'], bond[args.plot], label='Bond ' + bond['name'])
        
        for angle in angles:
            plt.plot(angle['x'], angle[args.plot], label='Bond ' + angle['name'])
            
        plt.legend(loc='upper right')
        plt.show()


if __name__ == '__main__':
    main()