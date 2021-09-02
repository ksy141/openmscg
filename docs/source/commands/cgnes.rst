cgnes
=====

.. automodule:: mscg.cli.cgnes

Notes
-----

* The command doesn't check the matrix sigularity for the normal equation. If the matrix contains a column of all zero values, the resulting coefficient for that column will be zero as well.
* Read `Ridge Regression (regularization) <https://en.wikipedia.org/wiki/Ridge_regression>`_ for more information about the Ridge model.
* Read the article `A Bayesian statistics approach to multiscale coarse graining <https://aip.scitation.org/doi/pdf/10.1063/1.3033218>`_ for more information about this approach.

Examples
--------

::
    
    cgnes --equation result/*.p --save combined.p
