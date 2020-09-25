Basic Concepts
==============

Introduction for basic concepts used in OpenMSCG.


Units
-----

In OpenMSCG, the following units are used:

* distance = Angstroms
* angle = Degree
* force = Kcal/mole/Angstrom
* temperature = Kelvin
* mass = grams/mole


Topology
--------

A topology file defines a group of particles (atoms or CG sites) and their bonding information for the simulated system. Different MD software packages have their own file format to define the topologies. For examples:

* LAMMPS - `LAMMPS Data Format <https://lammps.sandia.gov/doc/2001/data_format.html>`_

* GROMACS - `Toplogy file <http://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>`_

* CHARMM & NAMD - `Protein Structure File (PSF) <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`_

* AMBER - `PRMTOP File <https://ambermd.org/FileFormats.php#topology>`_

For the OpenMSCG package, we also have a self-customized file format to define the topologies of CG systems, named as **CGTOP**. Please read the detailed description of `CGTOP file format <cgtop.html>`_.


CLI Option for Topology
"""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a topology file for processing. A unified pattern for the option is as ``--top file``.

* The file format is automatically determined based on the content of the file.

* At this moment, only two different formats are supported by OpenMSCG, ``lammps`` and ``cgtop``. But more formats will be supported in future.

* In LAMMPS, the atom types are named in digital IDs (1, 2, ...). To be more user-friendly, the option ``--names`` can be used to reset the list of atom types by character-strings based names.

**Examples**::
    
    --top system.data
    --top system.top


Trajectory
----------

MD trajectories are collection of particle coordindates (also called as "frames") dumped from simulations. There are also many formats of trajectories from different MD software packages, such as `XYZ <https://en.wikipedia.org/wiki/XYZ_file_format>`_, `DCD <https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html>`_, `TRR <http://manual.gromacs.org/archive/5.0.3/online/xtc.html>`_ ...


CLI Option for Trajectory
"""""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a trajectory file for processing. A unified pattern for the option is as ``--traj file,[skip=n],[every=n],[frames=]``. The option value is consist of four fields separated by commas.

* ``file``: the path and name of the trajectory file.

The following optional fields are used for processing frames from a trajectory in a loop:

* ``skip``: skip over the first n frames for processing defined by ``skip``.
* ``every``: read a frame for every n frames defined by ``every``.
* ``frames``: read until n frames are read defined by ``frames``.

**Examples**::
    
    --traj md.trr,skip=10,every=100,frames=50

In this example, the file ``md.trr`` will be processed, while the first 10 frames will be skipped, then read a frame for every 100 frames (read one frame and then skip then next 99 frames), and the process loop will end after reading in a total of 50 frames.


Model and Table
---------------

In `molecular modelings <https://en.wikipedia.org/wiki/Force_field_(chemistry)>`_, models are mathematical functions of geometrical variables, i.e., distances or angles, to describe the forces and potential energies of a molecular system. One of the important steps in MSCG is to obtain the optimized parameters for the functions. Four styles of interactions can be specified as runtime options for either FM or REM methods: **--pair, --bond, --angle, --dihedral**.

CLI Option for Model
""""""""""""""""""""

For example, to fit model parameters for a pair-wise interaction between atom types A and B, one can specify the following option::

    --pair model=BSpline,type=A:B,min=3.0,max=8.0,resolution=0.2

The rumtime argument is followed by multiple ``key-value`` attributes separated by commas. Two attributes are mandatory:

* **model**: define the functional form. Currently only `BSpline` is allowed, but more functional types will be supported in future.
* **type**: a type name that used to specify the targeted interaction type to be optimized. In the example above, `A:B` indicates the pair-type that is formed by two atoms with types `A` and `B`.

Other attributes are optional and depend on the functional form. In the example above, there are three attributes to define a BSpline: **min, max and resolution**.



