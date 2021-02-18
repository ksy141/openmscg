Basic Concepts
==============

This chapter provides a basic outline of the objects and ideas used by openMSCG.


Units
-----

openMSCG uses the following units:

* distance = Angstroms
* angle = Degrees (for 3-body angles and 4-body dihedral torsions)
* force = Kcal/mol/Angstrom
* temperature = Kelvin
* mass = grams/mol


Topology
--------

A topology file defines a group of particles (atoms or CG sites) and their molecular bonding information.  Different MD software packages have their own file format to define topologies. For example:

* LAMMPS - `LAMMPS Data Format <https://lammps.sandia.gov/doc/2001/data_format.html>`_

* GROMACS - `Topology file <http://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>`_

* CHARMM & NAMD - `Protein Structure File (PSF) <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`_

* AMBER - `PRMTOP File <https://ambermd.org/FileFormats.php#topology>`_

For the openMSCG package, we offer a package-specific file format to define the topologies of CG systems, named **CGTOP**. The format description can be found `here <cgtop.html>`_.


Command Line Interface (CLI) Option for Topology
""""""""""""""""""""""""""""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a topology file for processing, and use the option flag ``--top file``.

* The file format (e.g. LAMMPS or CGTOP) is automatically determined based on the content of the file.

* At this moment, only two formats are supported by openMSCG, ``lammps`` and ``cgtop``. But more formats will be supported in the future.

* In LAMMPS, the atom types are named using numeric IDs (1, 2, ...). To be more user-friendly, the option ``--names`` can be used to define the list of atom types by character-string-based names.

**Examples**::
    
    --top system.data
    --top system.top


Trajectory
----------

MD trajectories are collections of particle coordinates (also called "frames") dumped from simulations. There are also many formats of trajectories from different MD software packages, such as `XYZ <https://en.wikipedia.org/wiki/XYZ_file_format>`_, `DCD <https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html>`_, `TRR <http://manual.gromacs.org/archive/5.0.3/online/xtc.html>`_.


CLI Option for Trajectory
"""""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a trajectory file for processing, and use the option flag ``--traj file,[skip=n],[every=n],[frames=]``. The option values consist of four fields separated by commas.

* ``file``: the path and name of the trajectory file.

The following optional fields are used for processing frames from a trajectory in a loop:

* ``skip``: skip over the first n frames for processing defined by ``skip``.
* ``every``: read a frame for every n frames defined by ``every``.
* ``frames``: read until n frames are read defined by ``frames``.

**Examples**::
    
    --traj md.trr,skip=10,every=100,frames=50

In this example, the file ``md.trr`` will be processed. First, the first 10 frames will be skipped, then 1 frame for every 100 frames (i.e. read 1 frame and then skip then next 99 frames) will be read, and the process will end after 50 frames have been read.


Models and Functional Forms
---------------------------

In molecular modeling, models (or force-fields) are mathematical functions of variables (e.g. distances or angles) to describe the forces and potential energies of a molecular system. One of the important outcomes in openMSCG is to obtain optimized parameters for the force-field. Four styles of interactions can be specified as runtime options for either FM or REM methods: **--pair, --bond, --angle, --dihedral**.



CLI Option for Model
""""""""""""""""""""

For example, to fit model parameters for a pairwise interaction using B-splines between atom types A and B, one can specify the following option::

    --pair model=BSpline,type=A:B,min=3.0,max=8.0,resolution=0.2

The runtime argument is followed by multiple ``key-value`` attributes separated by commas. Two attributes are mandatory:

* **model**: define the functional form. `BSpline` is the highlighted choice for multi-scale modeling, but new functional types are easy to be extended in **OpenMSCG**.
* **type**: name that is used to specify the targeted interaction types to be optimized. In the example above, `A:B` indicates a pair interaction between atoms of type `A` and `B`.

Other attributes are optional and depend on the functional form. In the example above, there are three attributes to define a BSpline: **min, max and resolution**.


List of Supported Models
""""""""""""""""""""""""

+----------+------------+--------------+---------------+
| Style    | Function   | FM supported | REM supported |
+==========+============+==============+===============+
| Pair     | BSpline    | Yes          | Yes           |
|          +------------+--------------+---------------+
|          | GaussCut   | No           | Yes           |
+----------+------------+--------------+---------------+
| Bond     | BSpline    | Yes          | Yes           |
+----------+------------+--------------+---------------+
| Angle    | BSpline    | Yes          | Yes           |
+----------+------------+--------------+---------------+
| Dihedral | BSpline    | Yes          | Yes           |
+----------+------------+--------------+---------------+


Force Field (Molecular Mechanics)
---------------------------------

Some of the functions/modules in OpenMSCG require calculating potential energies or forces for the system, which is defined by the conventional `force-field <https://en.wikipedia.org/wiki/Force_field_(chemistry)>`_ for pairs, bonds, angles ... OpenMSCG provides the ``force`` module for this calculation, which is usually specified by the option::

    --force force.yaml

where the parameter points to a `YAML` file as the definition of the force field. A sample file is as below::

    Pair_Table:
      CH3-CH3: Pair_CH3-CH3.table
      OH-OH:   Pair_OH-OH.table
      CH3-OH:  Pair_CH3-OH.table
    Bond_Harmonic:
      CH3-OH:  [5.0, 2.6]

The file contains a dictionary with multiple entries. The key of ach entry defines a type of potential function, and the value is another dictionary in which the keys are the targeted pair/bond/angle... types (as defined in the topology) and the values are the associated parameters.

The sample above demostrates the force field for a two-site methanol model including three types of pairwise interactions, using tabulated potential, and one type of bonded interaction, using the harmonic potential.

Supported Functional Forms
"""""""""""""""""""""""""""

+-------------------+---------------------------------------------------+
| Style             | Parameters                                        |
+===================+===================================================+
| Pair_Table        | Name of the table file                            |
+-------------------+---------------------------------------------------+
| Bond_Table        | Name of the table file                            |
+-------------------+---------------------------------------------------+
| Angle_Table       | Name of the table file                            |
+-------------------+---------------------------------------------------+
| Dihedral_Table    | Name of the table file                            |
+-------------------+---------------------------------------------------+
| Bond_Harmonic     | force constant, equilibrium value of the bond     |
+-------------------+---------------------------------------------------+
| Angle_Harmonic    | force constant, equilibrium value of the angele   |
+-------------------+---------------------------------------------------+
| Dihedral_Harmonic | force constant, equilibrium value of the dihedral |
+-------------------+---------------------------------------------------+

* The names of potential types are usually with two parts:

  * Type of interations: ``Pair``, ``Bond``, ``Angle`` ...
  
  * Type of funcional form: ``Table``, ``Harmonic`` ...
  
* For the table files, `LAMMPS format <https://lammps.sandia.gov/doc/pair_table.html>`_ is used.



