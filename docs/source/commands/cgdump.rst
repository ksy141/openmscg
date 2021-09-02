cgdump
======

.. automodule:: mscg.cli.cgdump

Notes
-----

This commands only supports LAMMPS tabulated potential at the moment, but support for other MD software packages will added in future.

The command generates tables with following attributes in the option:

* ``table_name``: the name of tabulated potential, in format of `style_type-type[-type]`. Examples: Pair_MeOH-MeOH, Bond_A-B, Angle_C1-C1-ET.
* ``xmin``: lower bound of output table.
* ``xmax``: upper bound of output table.
* ``padding`` (optional): If applying padding on the boundaries of table. The valaue can be `L` or `L2` - padding for lower boundary, `H` - padding for higher boundary, or `LH`, `L2H` - padding for both boundaries.

Padding rules:

* For lower bound, within the cut-offs in the **model**:
  
  * `L`: find the smallest i-th force value, where F\ :sub:`i` >0 and F\ :sub:`i` - F\ :sub:`i+1` >0, then use the gradient from (X\ :sub:`i`, F\ :sub:`i`) and (X\ :sub:`i+1`, F\ :sub:`i+1`) to expand the table from X\ :sub:`i` to X\ :sub:`xmin` by linear extrapolation.
  
  * `L2`: find the 3-rd and 4-th force values that are from the lower bound of the spline model, and calculate the slope from the two points. Then, use the slope to extraploate all other points lower than them.

* For upper bound, within the cut-offs in the **model**, find the largest i-th force value, where F\ :sub:`i` <0 and F\ :sub:`i-1` - F\ :sub:`i` >0, then use the gradient from (X\ :sub:`i-1`, F\ :sub:`i-1`) and (X\ :sub:`i`, F\ :sub:`i`) to expand the table from X\ :sub:`i` to X\ :sub:`xmax` by linear extrapolation.


Examples
--------

::
    
    cgdump --file result.p --dump Pair_MeOH-MeOH,2.5,8.0,0.1,

::
    
    cgdump --file result.p --plot Pair_MeOH-MeOH,2.5,8.0,0.1