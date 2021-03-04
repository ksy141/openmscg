cgib
====

.. include:: common.rst

.. automodule:: mscg.cli.cgib

* |cli-opt-top-name|

* |cli-opt-traj|

* The script generates a list of result files for each targeted variable specified in the options. The files are named as ``[Pair/Bond/Angle/Dihedral]-[Types].dat``. The result file contains three columns: the variable, the relative probability density, and the IB potential.

* The options ``--pair``, ``--bond``, and ``--angle`` are used to do histogram analysis for centain pair-types, bond-types, and angle-types. The multi-field values always comes with a number of type names (two for pair and bond, and three for angle), and are followed by the arguments for the histogram: min-value, max-value, and number of bins.

* If the option ``--plot`` is specified, the script will visualize the result by using the Python ``matplotlib`` package.


Notes
-----

* For the distribution of pair distances, the results have been normalized by the volume-correction factor, :math:`(4 \pi r^2)^{-1}`, and by average particle density probability so that the curve approach 1. These results corespond to the `radial distribution function (RDF) <https://en.wikipedia.org/wiki/Radial_distribution_function>`_.

* For other distributions (bonds, angles and dihedral torsions), results are scaled to give a maximum value of 1. So, the numbers can be intepreted as the relative probability to the most probable value.

* If you have used the old ``MSCGFM`` software, this command is a replacement to the ``rangerfinder``.


Examples
--------

::
    
    cgib --top tests/data/methanol_1728.data \
         --traj tests/data/methanol_1728_cg.trr \
         --names MeOH --cut 10.0 --temp 298.15 \
         --pair MeOH,MeOH,min=2.4,max=10.0,bins=200

::
    
    cgib --top tests/data/methanol_1728_2s.data \
         --traj tests/data/methanol_1728_2s.trr \
         --names CH3,OH --cut 10.0 --temp 298.15 \
         --pair CH3,CH3,min=2.4,max=10.0,bins=200 \
         --pair CH3,OH,min=2.4,max=10.0,bins=200 \
         --pair OH,OH,min=2.0,max=10.0,bins=200 \
         --bond CH3,OH,min=1.35,max=1.65,bins=60 \
         --plot n

    
