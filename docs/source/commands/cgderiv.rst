cgderiv
=======

.. include:: common.rst

.. automodule:: mscg.cli.cgderiv

* |cli-opt-top-name|

* |cli-opt-traj|

* |cli-opt-model|

Notes
-----

* The output of this script will be stored in a Python `Pickle file <https://docs.python.org/3/library/pickle.html#module-pickle>`_, which can be opened and read by standard Python envrionment. The data inside the Pickle file is packed as a dictionary with three attributes:

  * **models**: serialized model definitions in this task.
  * **dudl_mean**: ensemble averages of the derivatitives from each model stored as a dictionary indexed by the model name. The value of each model is a `Numpy` array.
  * **dudl_var**: similar to **dudl_mean**, but the contents are the variances of the derivatives.

Examples
--------

*To be added*

