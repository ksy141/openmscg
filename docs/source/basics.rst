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


Topology Files
--------------

A topology file defines a group of particles (atoms or CG sites) and their bonding information for the simulated system. Different MD software packages have their own file format to define the topologies. For examples:

* LAMMPS - `LAMMPS Data Format <https://lammps.sandia.gov/doc/2001/data_format.html>`_

* GROMACS - `Toplogy file <http://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html>`_

* CHARMM & NAMD - `Protein Structure File (PSF) <https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node23.html>`_

* AMBER - `PRMTOP File <https://ambermd.org/FileFormats.php#topology>`_

For the OpenMSCG package, we also have a self-customized file format to define the topologies of CG systems, named as **CGTOP**. Please read the detailed description of `CGTOP file format <cgtop.html>`_.


CLI Option for Topology
"""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a topology file for processing. A unified pattern for the option is as ``--top file-name,file-format``. The option value is consist of two fields separated by a comma.

* The first field is the path and name of the file, while the second field is the format of the file. The second field is required because it's hard to determine the format of the file automatically. (There's no universal standard for the file suffix name yet)

* At this moment, only two different formats are supported by OpenMSCG, ``lammps`` and ``cgtop``. But more formats will be supported in future.

* In LAMMPS, the atom types are named in digital IDs (1, 2, ...). To be more user-friendly, the option ``--names`` can be used to reset the list of atom types by character-strings based names.

**Examples**::
    
    --top system.data,lammps
    --top system.top,cgtop



Trajectory Files
----------------

MD trajectories are collection of particle coordindates (also called as "frames") dumped from simulations. There are also many formats of trajectories from different MD software packages, such as `XYZ <https://en.wikipedia.org/wiki/XYZ_file_format>`_, `DCD <https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html>`_, `TRR <http://manual.gromacs.org/archive/5.0.3/online/xtc.html>`_ ...


CLI Option for Trajectory
"""""""""""""""""""""""""

Many `CLI commands <commands.html>`_ need to read a trajectory file for processing. A unified pattern for the option is as ``--traj file-name,skip,every,count``. The option value is consist of four fields separated by commas.

* ``file-name``: the path and name of the trajectory file.

The following fields are used for processing frames from a trajectory in a loop:

* ``skip``: skip over the first n frames for processing defined by ``skip``.
* ``every``: read a frame for every n frames defined by ``every``.
* ``count``: read until n frames are read defined by ``count``.

**Examples**::
    
    --traj md.trr,10,100,50

In this example, the file ``md.trr`` will be processed, while the first 10 frames will be skipped, then read a frame for every 100 frames (read one frame and then skip then next 99 frames), and the process loop will end after reading in a total of 50 frames.


