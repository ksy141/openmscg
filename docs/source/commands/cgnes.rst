cgnes
=====

.. automodule:: mscg.cli.cgnes

Notes
-----

* Read `Ridge Regression (regularization) <https://en.wikipedia.org/wiki/Ridge_regression>`_ for more information.
* The command doesn't check the matrix sigularity for the normal equation. If the matrix contains a column of all zero values, the resulting coefficient for that column will be zero as well.

Examples
--------

::
    
    cgnes --equation result/*.p --save combined.p
