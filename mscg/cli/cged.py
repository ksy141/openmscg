'''CG Mapping for Essential Dynamics

Description
-----------

The ``cged`` is to optimize the mapping rule to capture the essential dynamics
from the all-atom trajectories (a.k.a., the EDCG method). The tool reads in the
result from a PCA analysis of the all-atom simulation trjectories, calculate
covariance matrix, and optmize the best grouping (mapping) rule, which gives
the minimal covariance residuals of the dynamics inside the grouped CG sites.

Usage
-----

Syntax of running ``cged`` command ::
    
    usage: cged.py [-h] [-v L] [--pc file] [--ev file] [--sites N] [--npc N]
                   [--save]

    Run EDCG optimization from PCA data. For detailed instructions please read
    https://software.rcc.uchicago.edu/mscg/docs/commands/cged.html

    General arguments:
      -h, --help         show this help message and exit
      -v L, --verbose L  screen verbose level (default: 0)

    Required arguments:
      --pc file          npy file for column-based 3xN eigenvectors (default:
                         None)
      --ev file          npy file for 3xN eigenvalues (default: None)
      --sites N          number of CG sites (N>1) (default: 2)

    Optional arguments:
      --npc N            number of principal components (default: 0)
      --save             file name for output (default: map)
      --cache            save residual matrix to file (default: False)
      
'''

from mscg import *
from time import time
import os

def compute_resmax(cov, with_cache=False):
    if os.path.isfile('resmax.chk.npy'):
        screen.info("Loading residual matrix from cache ...")
        return np.load('resmax.chk.npy').tolist()
    
    screen.info("Calculate residual matrix for all sub-domains ...")
    
    m = len(cov)
    res = [[0.0] * m for _ in range(m)]
    
    screen.info("Progress: ...", end="")
    count, last_count, time_start = 0, 0, time()
    
    
    # Algorithm with O(N^3)
    ''' 
    total = (m-1) * m * (2*m-1) / 12
    
    for sub_start in range(m-1):
        for sub_end in range(sub_start+1, m):
            res[sub_start][sub_end] = res[sub_start][sub_end - 1]
            L = list(range(sub_start, sub_end))
            
            for i in L:
                res[sub_start][sub_end] += cov[i][i] + cov[sub_end][sub_end] - 2.0 * cov[i][sub_end]
            
            count += len(L)
            res[sub_end][sub_start] = res[sub_start][sub_end]
        
            if count - last_count > 5000000:
                last_count = count

                elapsed = time() - time_start
                remained = (total - count) / count * elapsed
                screen.info("\rProgress: %0.2f%%. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (count*100/total, elapsed, remained), end="")
    '''
    
    # Algorithm with O(N^2)
    
    total = (m-1) * m / 2
    
    for L in range(1, m-1):
        for i in range(m-L):
            j = i + L
            
            res[i][j] = res[i][j-1] + res[i+1][j] - res[i+1][j-1]
            res[i][j] += cov[i][i] + cov[j][j] - 2.0 * cov[i][j]
            res[j][i] = res[i][j]
        
        count += m - L
        if count - last_count > 1000000:
            last_count = count

            elapsed = time() - time_start
            remained = (total - count) / count * elapsed
            screen.info("\rProgress: %0.2f%%. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (count*100/total, elapsed, remained), end="")
          
    screen.info("\nProgress: 100%%. Elapsed: %0.2f secs." % (time() - time_start))
    
    #print("\n".join([" ".join(["%10.4f" % (col) for col in row[:10]]) for row in res[:10]]))
    
    if with_cache:
        np.save('resmax.chk.npy',  np.array(res))
        
    return res

def minimize_dp(cov, n, with_cache=False):
    # residual matrix: res[i][j] - residuals for a sub-domain from atom i to j
    res = compute_resmax(cov, with_cache)
    
    # dp[i][j]: minimum residual for i splits for subdomain 0-j
    # dp[i][j] = min{dp[i-1][k-1] + cov[k][m], k=[i, m-1]}
    # best[i][j]: best position of the last split in i splits for subdomain 0-j
    # i = 0, n-1, the results for i splits (n domains <=> n-1 splits)
    
    m = len(res)
    dp = [[0.0] * m for _ in range(n)]
    best = [[None] * m for _ in range(n)]
    
    # dp[0][j] = no split, residual of subdomain from atom 0 to j
    dp[0] = res[0][:]
    
    # dp process for i split in domain [0 to j]
    total = sum([(m-i+1)*(m-i)/2 for i in range(1,n)])
        
    count, last_count, time_start = 0, 0, time()
    screen.info("Start global searching by dynamic programming (%d trials) ..." % (total))
    screen.info("Progress: ...", end="")
    
    for i in range(1, n):
        for j in range(i, m):
            # k -> split position: [0,k-1] and [k, j]
            positions = list(range(i, j+1))
            trials = [dp[i-1][k-1] + res[k][j] for k in positions]
            k = np.argmin(np.array(trials))
            dp[i][j] = trials[k]
            best[i][j] = k + i
            
            count += len(trials)
            
            if count - last_count > 2000000:
                last_count = count
                
                elapsed = time() - time_start
                remained = (total - count) / count * elapsed
                screen.info("\rProgress: %0.2f%%. Elapsed: %0.0f secs. Remaining %0.0f secs ..." % (count*100/total, elapsed, remained), end="")
    
    screen.info("\nProgress: 100%%. Elapsed: %0.2f secs." % (time() - time_start))
    
    # get recursive dp path
    split = [best[n-1][m-1]]
    
    for i in range(n-2, 0, -1):
        split.append(best[i][split[-1] - 1])
    
    return dp[n-1][m-1]/n/3, list(reversed(split))
    
    
    
def main(*args, **kwargs):
    
    # parse argument
    
    desc = 'Run EDCG optimization from PCA data. For detailed instructions please read ' + doc_root + 'commands/cged.html'
    
    parser = CLIParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter, fromfile_prefix_chars='@', add_help=False)
    
    group = parser.add_argument_group('General arguments')
    group.add_argument("-h", "--help", action="help", help="show this help message and exit")
    group.add_argument("-v", "--verbose", metavar='L', type=int, default=0, help="screen verbose level")
    
    group = parser.add_argument_group('Required arguments')
    group.add_argument("--pc", metavar='file', type=argparse.FileType('rb'), help="npy file for column-based 3xN eigenvectors")
    group.add_argument("--ev", metavar='file', type=argparse.FileType('rb'), help="npy file for 3xN eigenvalues")
    group.add_argument("--sites", metavar='N', type=int, default=2, help="number of CG sites (N>1)")
    
    group = parser.add_argument_group('Optional arguments')
    group.add_argument("--npc", metavar='N', type=int, default=0, help="number of principal components")
    group.add_argument("--save",  metavar='', type=str, default="map", help="file name for output")
    group.add_argument("--cache",  metavar='', type=bool, default=False, help="save residual matrix to file")
    
    if len(args)>0 or len(kwargs)>0:
        args = parser.parse_inline_args(*args, **kwargs)
    else:
        args = parser.parse_args()
    
    screen.verbose = args.verbose
    screen.info("OpenMSCG CLI Command: " + __name__)
        
    # read pca
    
    pc = np.load(args.pc)
    ev = np.load(args.ev)
    
    if len(pc.shape) != 2 or pc.shape[0] != pc.shape[1] or pc.shape[0]%3 != 0:
        raise Exception("Abnormal shape of the eigenvector matrix.")
    
    if pc.shape[0] != ev.shape[0]:
        raise Exception("Dimensions of eigenvectors and eigenvalues are not matched.")
    
    npc = ev.shape[0] if args.npc<=0 else (ev.shape[0] if args.npc>ev.shape[0] else args.npc)
    n = args.sites
    
    if n<2 or n>pc.shape[0] // 3:
        raise Exception("Incorrect number of CG sites.")
    
    # compute covariance matrix
    
    pc_dim = [pc[d::3, :npc].copy() for d in range(3)]
    ev = ev[:npc]
        
    cov_dim = [np.matmul(pc_dim[d], np.matmul(np.diag(ev), pc_dim[d].T)) for d in range(3)]
    cov = cov_dim[0] + cov_dim[1] + cov_dim[2]
    
    chi, splits = minimize_dp(cov.tolist(), n, args.cache)
    screen.info("Minimized Chi^2 (Sum of Residuals) = " + str(chi))
    
    # result
    
    segs_start = [0] + splits
    segs_end = [_-1 for _ in splits] + [pc.shape[0]//3-1]
    segs = np.array([segs_start, segs_end]).T
    
    if args.save == 'return':
        return chi, segs
    else:
        np.savetxt(args.save + ".txt", segs, fmt='%6.0f')
    
    
    
if __name__ == '__main__':
    main()
