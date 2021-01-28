'''Heterogeneous Elastic Network Models

Description
-----------

The ``cghenm`` is for building Heterogeneous Elastic Network Models (HENM).

Usage
-----

Syntax of running ``cghenm`` command ::

'''

from mscg import *
from mscg.sd import SD
from numpy import linalg as LA

def build_Hessian(N, bonds, K, r, v):
    H = np.zeros((N, N))
    
    for ib, bond in enumerate(bonds):
        a, b = bond[0], bond[1]
        for i in range(3):
            for j in range(3):
                sod = K[ib] * v[ib, i] * v[ib, j] / np.square(r[ib])
                H[a * 3 + i, b * 3 + j] = - sod
                H[b * 3 + j, a * 3 + i] = - sod
                
                H[a * 3 + i, a * 3 + j] += sod
                H[b * 3 + i, b * 3 + j] += sod
    
    return H
        

def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run normal equations solver for FM/UCG method. For detailed instructions please read ' + doc_root + 'commands/cgnes.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--mass", metavar='', type=float, default=[], nargs='*', help="masses of CG sites")
    group.add_argument("--temp", metavar='', type=float, default=300.0, help="temperature")
    group.add_argument("--alpha", metavar='', type=float, default=0.5, help="step size for iterations")
    group.add_argument("--maxiter", metavar='', type=int, default=1000, help="maximum iterations")
    group.add_argument("--ktol", metavar='', type=int, default=1.0e-4, help="tolerance of k-constants in iterations")
    
    group.add_argument("--sdstep", metavar='', type=float, default=1.0, help="initial step size for SD minimization")
    group.add_argument("--sdmax", metavar='', type=int, default=1000, help="maximum steps for SD minimization")
    group.add_argument("--sdftol", metavar='', type=float, default=1.0e-4, help="tolerance of forces SD minimization")
    
    group.add_argument("--traj", metavar='file[,args]', action=TrajReaderAction, help=TrajReaderAction.help, default=[])
        
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenCG CLI Command: " + __name__)
    
    # generate fluctuation info
    
    top, blist = None, None
    bonds = []
    
    for reader in TrajBatch(args.traj):
        
        if top is None:
            top = Topology()
            top.add_names('atom', ['CG'])
            top.add_names('bond', ['CG-CG'])
            top.add_atoms(['CG'] * reader.traj.natoms)
            
            ba1, ba2 = [], []
            
            for i in range(reader.traj.natoms):
                for j in range(i+1, reader.traj.natoms):
                    ba1.append(i)
                    ba2.append(j)
            
            top.add_bondings('bond', ['CG-CG'] * len(ba1), [ba1, ba2])
            blist = BondList(top.types_bond, top.bond_atoms, 
                             top.types_angle, top.angle_atoms, top.types_dihedral, top.dihedral_atoms)
        
        if reader.nread == 1:
            assert(reader.traj.natoms == top.n_atom)
        
        blist.build(reader.traj.box, reader.traj.x)
        bonds.append(blist.get_scalar('bond'))
        X, box = reader.traj.x.copy(), reader.traj.box.copy()
    
    bonds = np.stack(bonds)
    rf_mean = np.mean(bonds, axis=0)
    screen.info("MEAN(Bond-Length):" + str(rf_mean))
    rf_std = np.std(bonds, axis=0)
    screen.info("STDDEV(Bond-Length):" + str(rf_std))
    
    # iteration

    frc = Force({
        'Bond_Harmonic': {
            'CG-CG': [1.0, 1.0]
        }
    })
        
    frc.setup(top, None, blist)
    frc.funcs[0].v0 = rf_mean.copy()
    
    SDWorker = SD(args.sdstep, args.sdmax, args.sdftol)
    SDWorker.blist = blist
    SDWorker.force = frc
    
    N = X.shape[0]
    bonds = blist.bond
    v = np.zeros((len(bonds), 3), dtype=np.float32)
    K = np.ones(rf_mean.shape[0])
    
    for _ in range(args.maxiter):
        screen.info("\nStep %d ..." % _)
        
        # SD minimization
        frc.funcs[0].k = K.copy()
        X = SDWorker.run(X, box)
        
        # NMA
        blist.build(box, X)
        r = blist.get_scalar('bond')
        blist.get_vector('bond', v)
        
        # Build Hessian
        H = build_Hessian(N * 3, bonds, K * (2 * 4184 / args.temp / 1.38 / 6.022), r, v)
        
        # Eigen problem
        ew, ev = LA.eig(H)
        screen.info("Eigen-Values of Hessian Matrix: " + str(np.real(ew[:N*3-6])))
        
        # Vibration
        trial_std = np.zeros(len(bonds))
        
        for ib, bond in enumerate(bonds):
            a, b = bond[0], bond[1]
            prj = 0.0
            
            for imode in range(N * 3 - 6):
                sprj = 0.0
                
                for dim in range(3):
                    sprj += v[ib, dim] / r[ib] * (ev[a*3 + dim, imode] - ev[b*3 + dim, imode]) / np.sqrt(ew[imode])
            
                prj += np.square(sprj)
                    
            trial_std[ib] = np.sqrt(prj)
        
        screen.info("Current Fluctuation: " + str(trial_std))
        
        # Update K
        prevK = K.copy()
        K = 0.25 / (0.25 / K - args.alpha * (np.square(trial_std) - np.square(rf_std)))
        screen.info("New Constants: " + str(K))
        
        if np.abs(K - prevK).max() < args.ktol:
            break
    return
    
    

if __name__ == '__main__':
    main()
