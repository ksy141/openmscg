cgkmc
=====


*Available since 0.8.0*

Generate KMC-CG mappings of large biomolecules

Description
-----------

The `cgkmc` module realizes the implementation of K-Means Clustering Coarse-Graining (KMC-CG), a new method for generating optimal CG mappings for large biomolecules. KMC-CG removes the sequence-dependent constraints of ED-CG, allowing it to explore a more extensive space and thus enabling the discovery of more physically optimal CG mappings. Furthermore, the implementation of the K-means clustering algorithm can variationally optimize the CG mapping with efficiency and stability.


Usage
-----

Syntax of running `cgkmc` module

::

  Required arguments:
    beta        parameter controlling the positional residual
    gamma       parameter controlling the continuity residual
    frames      initialize the coordinates of CA atoms
  
  For more parameters, functions and how to run KMC-CG, please see the example file arp23.ipynb in the OpenMSCG/examples folder.


Notes
-----

The KMC-CG method reads in the aligned AA trajectories (after removing the transitional and rotational motions to the reference frame), and generates the optimal CG mappings based on the K-Means clustering method by minimizing a variational function that consists of a positional, a fluctuation, and a penalty term.

For more information about the KMC-CG method, see the following paper: `Wu, J.; Xue, W.; Voth, G. A. K-Means Clustering Coarse-Graining (KMC-CG): A Next Generation Methodology for Determining Optimal Coarse-Grained Mappings of Large Biomolecules. J. Chem. Theory Comput. 2023 <https://pubs.acs.org/doi/10.1021/acs.jctc.3c01053>`_.

A sample result is saved in the following format::
  
  0
  [0, 1, 2, 8, 9]
  1
  [3, 4, 5, 6, 7]

The number is the index of each CG site and the following list includes the residue indices belonging to that CG site. All the index starts from 0.