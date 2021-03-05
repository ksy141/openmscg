cghenm
======

*Available since 0.3.0*

.. include:: common.rst

.. automodule:: mscg.cli.cghenm

* |cli-opt-traj|


Notes
-----

The command reads in reference CG trajectories and calcualtes the mean and starndard deviation (referred as "flucuations") of the bonds.

From initial trial force-constants (``k=1.0``), the commands interate the following steps:

  1. Run the steepest descent algorithm to find the optimized structure (with the minimized potential energy)
  2. From the minimized structure, calculate trial fluctuations from the Hessian matrix
  3. Compare the trial fluctuations with the reference values and update the force constants ``k``.
  
A sample result is saved in the following format::

    Atom I     Atom J         R0          K    Fluc_Ref. Fluc_Matched
         1          2     15.938      9.540        0.180        0.180
         1          3     11.384     11.603        0.163        0.163
         2          3     23.743      1.521        0.450        0.450

The first two columns are the atoms of each bond, which is followed by the equlibrated bond lengths calculated as the mean values of bond lengths from the CG reference trajectories. The forth column contains the force constants fitted. The last two columns are the reference and matched bond fluctuations, respectively.

For more information about Heterogeneous elastic network models, see the following paper: `Lyman, E.; Pfaendtner, J.; Voth, G. A. Systematic Multiscale Parameterization of Heterogeneous Elastic Network Models of Proteins. Biophys. J. 2008, 95, 4183-4192 <http://www.sciencedirect.com/science/article/pii/S0006349508785588>`_.


Examples
--------

::
    
    cgfm --mass 5963 4355 3968 \
         --traj cg.lammpstrj \
         --temp 310.0 \
         --cut 30.0 \
         --alpha 0.5

