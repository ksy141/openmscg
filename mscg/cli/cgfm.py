
from mscg import *

def main(*args, **kwargs):
    
    # parse argument
    
    parser = CLIParser(description='Run MSCG force-matching method.', formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--top",  metavar='file,format', type=str, help="topology file: in format of [filename,format:lammps|cg]", required=True)
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--names",  metavar='', type=str, help="comma separated atom names (needed when using LAMMPS data file for topology)")
    group.add_argument("--traj", metavar='', type=str, help="trajectory files: in format of [filename,skip:0,every:1,nread:0]", action='append')
    group.add_argument("--cut", metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    group.add_argument("--pair", metavar='', type=str, help="add a pair RDF, in format of [type1,type2,min,res]", action='append')
    group.add_argument("--bond", metavar='', type=str, help="add a bond distribution, in format of [type1,type2,min,max,res]", action='append')
    group.add_argument("--angle", metavar='', type=str, help="add an angle distribution, in format of [type1,type2,type3,min,max,res]", action='append')
    group.add_argument("--save",  metavar='', type=str, default="matrix", help="file name for matrix output")
    
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
    
    # build up tables
    
    tables.empty()
    
    if args.pair is not None:
        for pair in args.pair:
            screen.info("Add pair coefficients table: " + pair)
            w = pair.split(',')
            TablePairBSpline(plist, '_'.join(w[0:2]), top.get_pair_type(w[0], w[1]), 
                             xmin=float(w[2]), resolution=float(w[3])).setup_cache()

    if args.bond is not None:
        for bond in args.bond:
            screen.info("Add bond coefficients table: " + bond)
            w = bond.split(',')
            TableBondBSpline(blist, '_'.join(w[0:2]), top.get_bond_type(w[0], w[1]), 
                             xmin=float(w[2]), xmax=float(w[3]), resolution=float(w[4])).setup_cache()

    if args.angle is not None:
        for angle in args.angle:
            screen.info("Add angle coefficients table: " + angle)
            w = angle.split(',')
            TableAngleBSpline(blist, '_'.join(w[0:3]), top.get_angle_type(w[0], w[1], w[2]), 
                              xmin=float(w[3]), xmax=float(w[4]), resolution=float(w[5])).setup_cache()
    
    # build up coefficients matrix
    
    screen.info("Build coefficients matrix ...")
    matrix = Matrix()
    matrix.add_tables(tables.all)
    matrix.setup(top.natoms)

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
            TIMER.click('io')
            TIMER.click('matrix', matrix.reset())
            TIMER.click('pair', plist.build(trj))
            TIMER.click('bond', blist.build(trj))
            TIMER.click('table', tables.compute_all())
            TIMER.click('matrix', matrix.multiplyadd(trj))

            # goto the next frame

            for i in range(every - 1):
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
        
    if args.save != "return":
        matrix.save("covariance_" + args.save)
    
    matrix.solve()
    
    if args.save != "return":
        matrix.save("coeffs_" + args.save)
    
    TIMER.click('solver')
    screen.info([""] + TIMER.report(False) + [""])
    
    # end
    
    if args.save == "return":
        return matrix.cov_y()
    

if __name__ == '__main__':
    main()
