''' Convert AA trajectories into a CG trajectory

Description
-----------

The `cgmap` command reads one or more all-atom (AA) simulation trajectories, and
converts the frames into a coarse-grained (CG) trajectory, according 
to the provided mapping rules.

Usage
-----

Syntax of running ``cgmap`` command ::

    usage: cgmap [-h] [-v] --map file --out file [--traj file[,args]]
    
    General arguments:
      -h, --help          show this help message and exit
      -v , --verbose      screen verbose level (default: 0)

    Required arguments:
      --map file          A YAML file with mapping rules (default: None)
      --out file          output trajectory (default: None)
      --traj file[,args]  reader for a trajectory file, multiple fields separated
                          by commas, the first field is the file name, while
                          others define the skip, every and frames (default args:
                          file,skip=0,every=1,frames=0) (default: [])
    
    Optional:
      --af, --no-af       Using the aggforce method for force-mapping (default: False)

Use the option "--af" will call the "aggforce" package that's developed by the AI4Science 
group to optimize (denoise) the CG forces with ML methods. When using "aggforce", there
are another set of arguments that can be set. Please use "cgmap -h" to gain the full list
of arguments, or visit the link below for the official documentation of this method:

* https://github.com/noegroup/aggforce
'''

from mscg import *

def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run MSCG Mapping for generation of CG trajectories. For detailed instructions please read ' + doc_root + 'commands/cgmap.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--map", metavar='file', type=str, required=True, help="A YAML file with mapping rules")
    group.add_argument("--out", metavar='file', type=str, help="output trajectory", required=True)
    group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])
    
    _aggforce_help = 'See aggforce documentation'

    group.add_argument("--af", type=bool, default = False, action=argparse.BooleanOptionalAction, help="Using the aggforce method for force-mapping")
    group.add_argument("--af-map-style", type=str, default='basic', help=_aggforce_help)
    group.add_argument("--af-constraint-threshold", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-constraint-frames", type=int, default=10, help=_aggforce_help)
    group.add_argument("--af-inner", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-outer", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-width", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-n-basis", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-batch-size", type=float, default=1e-3, help=_aggforce_help)
    group.add_argument("--af-lazy", type=bool, default=True, help=_aggforce_help)
    group.add_argument("--af-temperature", type=float, default=300, help=_aggforce_help)
    group.add_argument("--af-l2", type=float, default=1e-3, help=_aggforce_help)

    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)
    
    # building map
    
    mapper = Mapper.build_from_yaml(args.map)
    site_types = mapper.get_types()
    
    # AGGFORCE mapping

    if args.af:
        screen.info("The method aggforce mapping is used. Importing the package now ...")
        
        try:
            #import aggforce
            from aggforce import linearmap as lm
            from aggforce import agg as ag
            from aggforce import constfinder as cf
            from aggforce import featlinearmap as p
            from aggforce import jaxfeat as jf
        except:
            screen.fatal("failed to import the aggforce package. Please check the installation.")
        
        TIMER.reset()
        screen.info("Load all force data for denoising ...")
        x, f = [], []

        for reader in TrajBatch(args.traj):
            assert reader.traj.has_attr('f')
            f.append(reader.traj.f.copy())

            if  args.af_map_style == 'optim_config' or len(x) < args.af_constraint_frames:
                 x.append(reader.traj.x.copy())

            TIMER.click('agg_load')
        
        map_matrix = mapper.get_matrix(f[-1].shape[0])
        cmap = lm.LinearMap(map_matrix)

        screen.info("Finding constraints...")

        constraints = cf.guess_pairwise_constraints(
            np.array(x[0:args.af_constraint_frames]),  
            threshold = args.af_constraint_threshold)
        
        TIMER.click('agg_cnstr')
        screen.info("Optimizing force mapping...")
        
        def project_forces(method, **kwargs):
            return ag.project_forces(xyz = None, forces = np.array(f), 
                config_mapping = cmap, constrained_inds = constraints, method = method)
        
        if args.af_map_style == 'basic':
            optim_results = project_forces(lm.constraint_aware_uni_map)
        elif args.af_map_style == 'optim':
            optim_results = project_forces(lm.qp_linear_map)
        elif args.af_map_style == 'optim_config':
            config_feater = p.Curry(jf.gb_feat, inner = args.af_inner, outer = args.af_outer, width = args.af_width, 
                n_basis = args.af_n_basis, batch_size = args.af_batch_size, lazy = args.af_lazy)
            
            optim_results = ag.project_forces(xyz = np.array(x), forces = np.array(f), config_mapping = cmap, 
                constrained_inds = constraints, l2_regularization = 1e3, kbt = 0.001987204259 * args.af_temperature,
                featurizer = p.Multifeaturize([p.id_feat, config_feater]), method = p.qp_feat_linear_map)
            
        else:
            screen.fatal("Illegal map style: [%s]." % (args.af_map_style))

        F_optim = list(optim_results['projected_forces'])
        print(F_optim)

        TIMER.click('agg_optim')
        screen.info("Force mapping optimization is finished.")
        screen.info([""] + TIMER.report(False) + [""])
        
    # start processing trajectory
    
    if args.out == 'return':
        frames = {'box':[], 'type':site_types, 'x':[], 'f':[]}
    else:
        outfile = Trajectory(args.out, "w")
        outfile.t = np.array(site_types)
        outfile.v = None
    
    TIMER.reset()
    last = TIMER.last
    
    for reader in TrajBatch(args.traj):
        TIMER.click('read')
        
        X, F = mapper.process(reader.traj.box, reader.traj.x, reader.traj.f if reader.traj.has_attr('f') and args.af is False else None)
        TIMER.click('map')

        if args.af:
            F = F_optim.pop()
        
        if args.out == 'return':
            frames['box'].append(reader.traj.box.copy())
            frames['x'].append(X)
            frames['f'].append(F)
        else:
            outfile.box, outfile.x, outfile.f = reader.traj.box.copy(), X, F
            outfile.timestep = reader.traj.timestep
            outfile.write_frame()
        
        TIMER.click('write')
        
    screen.info("Processing is finished.")
    screen.info([""] + TIMER.report(False) + [""])
    
    if args.out == 'return':
        return frames


if __name__ == '__main__':
    main()
