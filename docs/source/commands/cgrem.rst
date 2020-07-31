cgrem
=====

This command is used to dump the force tables for different MD software packages from the result generated from the ``cgfm`` command.


Usage
-----

Syntax of running ``cgdump`` command ::

    usage: cgdump [-h] [-v L] --file file [--dump] [--plot]
    
    General arguments:
      -h, --help         show this help message and exit
      -v L, --verbose L  screen verbose level (default: 0)

    Required arguments:
      --file file        matrix file (default: None)

    Optional arguments:
      --dump             dump a LAMMPS table file, in format of
                         [table_name,xmin,xmax,dx] (default: None)
      --plot             plot table, in format of [table_name,xmin,xmax,dx]
                         (default: None)


Notes
-----

This commands supports only the LAMMPS tabulated potential at the moment, but supports to other MD software packages will added in future.


Examples
--------

*To be added*
