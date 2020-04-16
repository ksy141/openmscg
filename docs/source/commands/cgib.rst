cgib
====

Description
-----------

The ``cgib`` command is used to analyze the MD trajectoris and calculate the distribution histograms for targeted structural variables: **pairs, bonds, angles and dihedral torsions**. The function can be used to find the ranges of these variables, which can be used for the force-matching method later. This command also calculates the **Boltzmann-Inversed (IB)** free energy profiles for the given variables.


Usage
-----

Syntax of running ``cgib`` command ::

    usage: cgib [-h] [-v] --top file [--names] [--traj file[,args]] [--cut]
            [--temp] [--pair types,args] [--bond types,args]
            [--angle types,args] [--plot U or N]

    General arguments:
      -h, --help          show this help message and exit
      -v , --verbose      screen verbose level (default: 0)

    Required arguments:
      --top file          topology file (default: None)

    Optional arguments:
      --names             comma separated atom type names (needed when using
                          LAMMPS data file for topology) (default: None)
      --traj file[,args]  reader for a trajectory file, multiple fields separated
                          by commas, the first field is the file name, while
                          others define the skip, every and frames (default args:
                          file,skip=0,every=1,frames=0) (default: [])
      --cut               cut-off for pair interactions (default: 10.0)
      --temp              temperature (K) for IB (default: 298.15)
      --pair types,args   define new pair analysis with format: type1,type2,args;
                          args and default values are: min=0,max=10,bins=10
                          (default: [])
      --bond types,args   define new bond analysis with format: type1,type2,args;
                          args and default values are: min=0,max=10,bins=10
                          (default: [])
      --angle types,args  define new angle analysis with format:
                          type1,type2,type3,args; args and default values are:
                          min=0,max=10,bins=10 (default: [])
      --plot U or N       plot the results of U (potential) or n (distribition)
                          (default: U)



* For more details about the option ``--top`` and ``--name``, check the `topology concept <../basics.html#cli-option-for-topology>`_.

* For more details about the option ``--traj``, check the `trajectory concept <../basics.html#cli-option-for-trajectory>`_.

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