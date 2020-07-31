Commands
========

This document details the basic usages of OpenMSCG for the work of multi-scale modeling. These functions are warpped as command-line-interface (CLI) scripts installed in the ``bin`` folder of ``Python`` or ``Conda``. Users can directly execute these scripts in the command-line terminals.


Command-Line Options
--------------------

All the CLI scripts needs to be executed with one or more options (`General concepts about CLI options <https://docs.python.org/3/howto/argparse.html>`_). These options are used to pass the input parameters to the task. To get a full list of options, users can simply execute it with ``-h`` or ``--help``. For example ::

    $ cgfm -h
    usage: cgfm [-h] [-v L] --top file [--names] [--traj file[,args]] [--cut]
            [--save] [--pair types,args] [--bond types,args]
            [--angle types,args]

    Run MSCG force-matching method. For detailed instructions please read
    https://software.rcc.uchicago.edu/mscg/docs/commands/cgfm.html

    General arguments:
      -h, --help          show this help message and exit
      -v L, --verbose L   screen verbose level (default: 0)

    Required arguments:
      --top file          topology file (default: None)

    Optional arguments:
      --names             comma separated atom names (needed when using LAMMPS
                          data file for topology) (default: None)
    
    ...


* The values of an option is specified in the format of ``--option values``. Some of the options also have a shortened name after a single bar, i.e., ``--verbose`` can be also used as ``-v``.

* Some options are required, which must be passed to execute the script; while the others are **optional**. If the optional option is not specified, a default value will be used.

* For some options, the value has multiple fields, for example, the option ``--traj`` for the command ``cgfm`` takes one or more fields, where there is at least one field for the name of the trajectory file and may have other optional fields. These fields are separated by comma. But, no space is allowed in the whole character-string of the value.

.. admonition:: Note
  
  For an multi-field value, it is started with several fields as positional arguments (value only), and (or) followed by other ones as optional arguments (``key=value`` format). And, the optional arguments have default values to be used, if they are not specified in the command. For example, the option ``--traj`` has the format of value as:
  
  ``--traj file,skip=0,every=1,frames=0``
  
  where the positional field *file* is the name of the file, and the optional arguments are *skip* (number of frames skipped at the beginning), *every* (read one frame every this many frames) and *frames* (max number of frames will be read).
  

* Some options can be specified multiple times, resulting in a list of values passed for these options. For example, ``--bond`` can be specified multiple times to add more than one target bond tables to the force-matching algorithm.


Take Options from File
----------------------

Usually the scripts require many options, so users may need to explictly typing all options to form a really long command-line, which is very innconenient. Python provides an alternative to simplify it: **taking options from a text file**.

For example, user can put options in a file (i.e., ``cgib.opts``) like ::

    --top run/data/lmp.data,lammps
    --traj run/data/cg_methanol_1728_2s.trr
    --traj run/data/cg_methanol_1728_2s.trr,every=50
    --names CH3,OH
    --cut 10.0
    --temp 298.15
    --pair CH3,CH3,min=2.8,max=10.0,bins=200
    --pair CH3,OH,min=2.8,max=10.0,bins=200
    --pair OH,OH,min=2.5,max=10.0,bins=200
    --bond CH3,OH,min=1.35,max=1.65,bins=60
    --plot none
    --verbose 1

Then, the script can be executed by passing this file to it after the "@" symbol ::
    
    $ cgib @chib.opts

Also, the two ways of specifying options (from command-line or a file) can be mixed, which provide a practical solution to conduct a batch of tasks. For example, if we can put the options with the same values in an option file, while leaving the variable options in the command-line when executing it.


Call Command in Python
----------------------

Finally, these script can be also directly called from Python codes. Each script is actually a Python module in the ``cli`` sub-package of OpenMSCG. Executing the script in the CLI is equivalent as calling the ``main`` function in that module. For example, we can execute the script ``cgib`` like ::
    
    $ cgib --top run/data/lmp.data,lammps \
      --traj run/data/cg_methanol_1728_2s.trr \
      --pair CH3,OH,min=2.8,max=10.0,bins=200 \
      --pair OH,OH,min=2.5,max=10.0,bins=200 \
      --names CH3,OH

which is equivalent as following codes in Python::
    
    from mscg.cli import cgib
    
    cgib.main(top='run/data/lmp.data,lammps',
        traj='run/data/cg_methanol_1728_2s.trr',
        pair=['CH3,OH,min=2.8,max=10.0,bins=200', 'OH,OH,min=2.5,max=10.0,bins=200'],
        names='CH3,OH')

* These options can be interpreted as the `**kwargs <https://book.pythontips.com/en/latest/args_and_kwargs.html>`_ when calling a function.

* For the options that can be specified multiple times, i.e., ``--pair``, the values can be aggregated as a list passed to the argument of the function.


List of Commands
----------------

.. toctree::
   :maxdepth: 1
   
   commands/cgib
   commands/cgmap
   commands/cgfm
   commands/cgdump
   commands/cgderiv
   commands/cgrem
