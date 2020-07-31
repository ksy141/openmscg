
from mscg import *

def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run MSCG Mapping for generation of CG trajectories. For detailed instructions please read ' + doc_root + 'commands/cgmap.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--map", metavar='file', type=str, required=True)
    group.add_argument("--out", metavar='file', type=str, help="output trajectory", required=True)
    group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])
        
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenCG CLI Command: " + __name__)
    
    # building map
    
    mapper = Mapper.build_from_yaml(args.map)
    
    # start processing trajectory
    
    if args.out == 'return':
        frames = {'box':[], 'x':[], 'f':[]}
    else:
        outfile = Trajectory(args.out, "w")
        outfile.type = outfile.v = None
    
    TIMER.reset()
    last = TIMER.last
    
    for reader in TrajBatch(args.traj):
        TIMER.click('read')
        
        X, F = mapper.process(reader.traj.box, reader.traj.x, reader.traj.f if reader.traj.has_force() else None)
        TIMER.click('map')
        
        if args.out == 'return':
            frames['box'].append(reader.traj.box.copy())
            frames['x'].append(X)
            frames['f'].append(F)
        else:
            outfile.box, outfile.x, outfile.f = reader.traj.box.copy(), X, F
            outfile.write_frame()
            
        TIMER.click('write')
        
    screen.info("Processing is finished.")
    screen.info([""] + TIMER.report(False) + [""])
    
    if args.out == 'return':
        return frames


if __name__ == '__main__':
    main()