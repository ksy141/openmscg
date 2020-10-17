cgdump
======

.. automodule:: mscg.cli.cgdump

Notes
-----

This commands only supports LAMMPS tabulated potential at the moment, but support for other MD software packages will added in future.


Examples
--------

::
    
    cgdump --file result.p --dump Pair_MeOH-MeOH,2.5,8.0,0.1

::
    
    cgdump --file result.p --plot Pair_MeOH-MeOH,2.5,8.0,0.1