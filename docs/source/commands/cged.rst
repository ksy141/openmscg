cged
====

*Available since 0.3.2-dev*

.. automodule:: mscg.cli.cged

Notes
-----

Full details about the EDCG method is discussed in this `article <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2586547/>`_.

Principal Component Analysis
""""""""""""""""""""""""""""

* The mapping essential dynamics is always started from principal component analysis (PCA) on the all-atom simulation trajectories. However, this step is **NOT** included in this tool, because there are already lots of software packages that can do this step with high performances, such as:

    * `GROMACS <https://manual.gromacs.org/documentation/2019-rc1/reference-manual/analysis/dihedral-pca.html>`_
    * `VMD/NWWIZ <https://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/>`_
    * `MDTraj <https://mdtraj.org/1.9.4/examples/pca.html>`_
    * `MDAnalysis <https://docs.mdanalysis.org/1.0.0/documentation_pages/analysis/pca.html>`_

  An important prior step before PCA is that the translation and angular 
  rotations of the full system needs to be eliminated from the trajectories.
  This step is usually referred as to `alighment` or `fit-to` functions in
  the software tools. An example of Python script of PCA on the C-alpha atoms
  is below::
      
      import MDAnalysis as mda
      from MDAnalysis.analysis import align, pca
      
      u_traj = mda.Universe('protein.psf', 'protein.dcd')
      u_ref = mda.Universe('protein.psf', 'protein.dcd')

      u_traj.trajectory[-1]
      u_ref.trajectory[0]
      
      align.AlignTraj(u_traj, u_ref, in_memory=True).run()
      pc = pca.PCA(u_traj, align=False, mean=None, n_components=None).run()
      
      import numpy as np
      np.save('pc.npy', pc.p_components)
      np.save('ev.npy', pc.variance)


* The tool reads in two data sets from the PCA results, eigenvectors and eigenvalues, both of which should be prepared in `NumPy` binary format (npy). Eigenvalues are stored in a 1-d NumPy array sorted in descending order, while eigenvectors are stored in a 2-d NumPy array, and each column is an eigenvector that matched the order of eigenvalues.

* After reading the PCA data, the covariance matrix will be calculated. Theoption ``--npc`` defines the number of the largest principal components will be used in this step (a value 0 means using all components).


Minimization Algorithm and Performance
""""""""""""""""""""""""""""""""""""""

In the original article the minimized covariance residual can be found by 
**Simulated Annealing (SA)** and **Steepest Descent (SD)**, but none of these
methods can ensure the global minimum. Therefore, in OpenMSCG we introduced the
**Dynamic Programming** (`DP <https://en.wikipedia.org/wiki/Dynamic_programming>`_)
algorithm implementation, which is able to find the global minimum like the
Depth-First-Search (DFS), but much more efficient than it.

DP requires calculating the residual matrix, ``resmax[i][j]``, which stores
the residuals for all sub-segments from i-th atom to j-th atom. The time-complexity
of calculating this matrix is O(N^2), where N is the number of atoms. After 
calculation of residual matrix, the time-complexity for searching global minimum i
s O(N^2 x P), where P is the number groups (CG sites) to be determined.

In real practices, people may try different numbers of CG sites from the same
covariance matrix. Therefore, a cache file named as ``resmax.cache.npy`` is
created in the same folder storing the data of the residual matrix. When rerunning
the tool with different ``--sites`` option, the residual matrix can be read from
this cache file instead of computing it again. However, if the PCA data are
updated, or the option ``--npc`` is changed, users need to manually remove and
re-generate this cache file.

.. admonition:: Note

  In a testing of CGED on `Intel(R) Xeon(R) Gold 6148 CPU @ 2.40GHz` to map a system 
  with N=372 alpha carbons into P=10 CG sites:
  
  * The time of calculating residual matrix is about `0.07` seconds, which means it requires about 50 seconds for a system of 10,000 alpha carbons.
  
  * The time of optimization takes about `0.33` seconds for 10 CG sites, which means it requires about 4 minutes for a system of 10,000 alpha carbons. And, mapping this system to 1000 CG sites will require about 400 minutes (~7 hours).
  
  Therefore, the time of this algorithm is totally acceptable for systems with less
  than 10k sites to be mapped.

The DP altorighm implemented currently has two limits:

1. it can be used only for sequency-based mapping (consecutive grouping).
2. not suitable for mapping of more than 10k atoms.

Multithreading accerlation and SA/SD algorithms will be developed in future versions.



Examples
--------

::
    
    cged --pc pc.npy \
         --ev ev.npy \
         --npc 24 \
         --sites 10

    