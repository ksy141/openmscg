cgib
====

.. include:: common.rst

.. automodule:: mscg.cli.cgib

* |cli-opt-top-name|

* |cli-opt-traj|

* The script generates a list of result files for each targeted variable specified in the options. The files are named as ``[Pair/Bond/Angle/Dihedral]-[Types].dat``. The result file contains three columns as the variable, relative probability density and IB potential.

* The options ``--pair``, ``--bond`` and ``--angle`` are used to do a histogram analysis for centain pair-types, bond-types and angle-types. The multi-field values are always comes with a number of type names (two for pair and bond, and three for angle), and followed by the arguments for the histogram: min-value, max-value and number of bins.

* If the option ``--plot`` is specified, the script will visualize the result by using the Python ``matplotlib`` package.


Notes
-----

* For the distribution of pair distances, the results are already applied with the volume-correction factor, :math:`(4 \pi r^2)^{-1}`, and to the average particle density probability as 1. Therefore, it is the comonly called `radial distribution function (RDF) <https://en.wikipedia.org/wiki/Radial_distribution_function>`_.

* For other distributions (bonds, angles and dihedral torsions), the results are just simply scaled to give the maximum value of 1. So, the numbers can be intepreted as the relative probability to the most probable value.

* If you have used the old ``MSCGFM`` software, this command is a replacement to the ``rangerfinder`` command in ``MSCGFM``.


Examples
--------

*To be added*
