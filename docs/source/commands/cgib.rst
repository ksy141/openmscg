cgib
====

The `cgib` command performs non-iterative Boltzmann Inversion.

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

*To be added*
