'''Force-matching to construct force tables with B-Spline or other linear models

Description
-----------

The ``cgmd`` is a mini simulator.

Usage
-----

Syntax of running ``cgfm`` command ::

    usage: cgfm.py [-h] [-v L] --top file [--names] [--traj file[,args]] [--cut]
                   [--save] [--alpha] [--bayesian] [--pair [key=value]]
                   [--bond [key=value]] [--angle [key=value]]
                   [--dihedral [key=value]] [--ucg [key=value]]
                   [--ucg-wf [key=value]] [--force]

    General arguments:
      -h, --help            show this help message and exit
      -v L, --verbose L     screen verbose level (default: 0)

    Required arguments:
      --top file            topology file (default: None)

    Optional arguments:

'''

from mscg import *

def main(*args, **kwargs):

    # parse argument

    desc = 'Run CG simulation. For detailed instructions please read ' + doc_root + 'commands/cgmd.html'

    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)

    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")

    group = parser.add_argument_group('Required arguments')
    group.add_argument("--top",  metavar='file', action=TopAction, help="topology file", required=True)
    group.add_argument("--start",  metavar='file', action=TrajReaderAction, help="starting coordinates", default=[], required=True)

    group.add_argument("--force",  metavar='file', type=argparse.FileType('r'), help="force-field definitions", required=True)
    group.add_argument("--cut", metavar='', type=float, default=10.0, help="cut-off for pair interactions")
    group.add_argument("--exclude", metavar='', type=str, default="111", help="exclude 1-2, 1-3, 1-4 bonding neighbors in the pair-list")

    #group = parser.add_argument_group('Optional arguments')
    #

    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()

    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)

    # Read topology information

    screen.info("Read topology ... ")

    #if args.names is not None:
    #    args.top.reset_names(args.names.split(','))

    # Prepare lists

    screen.info("Build pair and bonding lists ...")
    plist = PairList(cut = args.cut, binsize = args.cut * 0.5)
    plist.init(args.top.types_atom, args.top.linking_map(*([bit=='1' for bit in args.exclude[:3]])))
    blist = BondList(
        args.top.types_bond, args.top.bond_atoms,
        args.top.types_angle, args.top.angle_atoms,
        args.top.types_dihedral, args.top.dihedral_atoms)

    # Read force-field information

    screen.info("Read force-field ...")
    force = Force(yaml.load(args.force, Loader=yaml.FullLoader))
    force.setup(args.top, plist, blist)

    # Force

    args.start[0].next_frame()
    start = args.start[0].traj

    plist.setup_bins(start.box)
    plist.build(start.x)
    blist.build(args.start[0].traj.box, args.start[0].traj.x)
    force.compute()
    print(force.f)

if __name__ == '__main__':
    main()
