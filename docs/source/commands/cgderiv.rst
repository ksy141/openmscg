cgderiv
=======

.. include:: common.rst

.. automodule:: mscg.cli.cgderiv

* |cli-opt-top-name|

* |cli-opt-traj|

* |cli-opt-model|

Notes
-----

* The output of this script will be stored in a Python `Pickle file <https://docs.python.org/3/library/pickle.html#module-pickle>`_, which can be opened and read by standard Python libraries. The data inside the Pickle file is a dictionary with three keys:

  * **models**: serialized model definitions.
  * **dudl_mean**: ensemble averages of the potential energy parameterization derivatives from each model stored as a dictionary indexed by the model name. The value of corresponding to each model is a `Numpy` array.
  * **dudl_var**: similar to **dudl_mean**, but the contents are the variances of the derivatives.

Examples
--------

*To be added*

