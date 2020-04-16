cgfm
====


The ``cgib`` command is used to do the force-matching method with a given CG topology and a set of CG trajectories. It is the main computing engine of the OpenMSCG package. Briefly, it reads in the coordinates and forces from the trajectories, calculating the loadings for a group of user-defined force tables, storing them in a covariance matrix, and finaly solving the linear regression to get the coefficients for each force types. Based on the variational theory, the computed coefficients defines the force function that minimize the force residues between real CG values and models. 

In OpenMSCG, a force table is defined in B-spline function defined by *n* **uniformly** spaced knots. Users can specify the starting value (``min``), ending value (``max``), the interval between adjacent knonts (``resolution``) and the ``order`` of the spline. Finally, the number of knots, which is also the number of computed coefficients, is

``n = (max - min) / resolution + order - 2``


Usage
-----

Syntax of running ``cgib`` command ::

    usage: cgfm [-h] [-v L] --top file [--names] [--traj file[,args]] [--cut]
                [--save] [--pair types,args] [--bond types,args]
                [--angle types,args]
    
    General arguments:
      -h, --help          show this help message and exit
      -v L, --verbose L   screen verbose level (default: 0)

    Required arguments:
      --top file          topology file (default: None)

    Optional arguments:
      --names             comma separated atom names (needed when using LAMMPS
                          data file for topology) (default: None)
      --traj file[,args]  reader for a trajectory file, multiple fields separated
                          by commas, the first field is the file name, while
                          others define the skip, every and frames (default args:
                          file,skip=0,every=1,frames=0) (default: [])
      --cut               cut-off for pair interactions (default: 10.0)
      --save              file name for matrix output (default: matrix)
      --pair types,args   define new pair table with format: type1,type2,args:
                          min,max,resolution,order (default: [])
      --bond types,args   define new bond table with format: type1,type2,args:
                          min,max,resolution,order (default: [])
      --angle types,args  define new angle table with format:
                          type1,type2,type3,args: min,max,resolution,order
                          (default: [])


* For more details about the option ``--top`` and ``--name``, check the `topology concept <../basics.html#cli-option-for-topology>`_.

* For more details about the option ``--traj``, check the `trajectory concept <../basics.html#cli-option-for-trajectory>`_.

* Users can use options ``--pair``, ``--bond``, ``--angle`` to add force tables to be computed. These options are followed by multi-field values, which comes first with a number of atom types defining the type of forces and then several optional arguments defining the B-Spline knots as described above.


Notes
-----

* The script ``cgfm`` doesn't dump the force or coefficient table in the plain-text format. Instead, there are two binary output files, in the `Python pickle <https://docs.python.org/3/library/pickle.html>`_ format, named by the option ``--save``. By default, they are *covariance_matrix.p* and *coeffs_matrix.p*.

  * **coeffs_matrix.p** stores the coeffients of B-Spline for each force table after solving the regression problem. Users can use the command `cgdump <cgdump.html>`_ to dump the readable force tables from it.
  
  * **covariance_matrix.p** stores the covariance matrix for the regression problem before solving it. It can be used as intermediate results for debugging or other purposes.

* If this command is `called with in another Python script <../commands.html#call-command-in-python>`_, users may want to access the results after call its main function. In this case, users can pass the value *"return"* to the ``--save`` function, and the coefficient matrix will be returned to the wrapper script as lists. Example::
    
    coeffs = cgfm.main(
        top     = datafile("unary_lj_fluid.top"),
        traj    = datafile("unary_lj_fluid.lammpstrj,frames=20"),
        cut     = 2.50,
        pair    = ['1,1,min=0.9,resolution=0.1'],
        save    = 'return'
    )

* If you have used the old ``MSCGFM`` software, this command is a replacement to the ``newfm`` command in ``MSCGFM``.


Examples
--------

*To be added*