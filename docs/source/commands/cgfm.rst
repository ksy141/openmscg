cgfm
====

.. include:: common.rst

.. automodule:: mscg.cli.cgfm


* |cli-opt-top-name|

* |cli-opt-traj|

* Users can use options ``--pair``, ``--bond``, ``--angle`` to add force tables to be computed. These options are followed by multi-field values, which comes first with a number of atom types defining the type of forces and then several optional arguments defining the B-Spline knots as described above.

* For more details about the options ``--ucg`` and ``--ucg-wf``, please read the section for `UCG Method`_.

Notes
-----

* The script ``cgfm`` doesn't dump the force or coefficient table in the plain-text format. Instead, there are two binary output files, in the `Python pickle <https://docs.python.org/3/library/pickle.html>`_ format, named by the option ``--save``. By default, they are *covariance_matrix.p* and *coeffs_matrix.p*.

  * **coeffs_matrix.p** stores the coeffients of B-Spline for each force table after solving the regression problem. Users can use the command `cgdump <cgdump.html>`_ to dump the readable force tables from it.
  
  * **covariance_matrix.p** stores the covariance matrix for the regression problem before solving it. It can be used as intermediate results for debugging or other purposes.

* If this command is `called with in another Python script <../commands.html#call-command-in-python>`_, users may want to access the results after call its main function. In this case, users can pass the value *"return"* to the ``--save`` function, and the coefficient matrix will be returned to the wrapper script as lists. Example::
    
    coeffs = cgfm.main(
        top     = datafile("unary_lj_fluid.top"),
        traj    = datafile("unary_lj_fluid.lammpstrj,frames=20"),
        cut     = 2.50,
        pair    = ['1,1,min=0.9,resolution=0.1'],
        save    = 'return'
    )

* If you have used the old ``MSCGFM`` software, this command is a replacement to the ``newfm`` command in ``MSCGFM``.


UCG Method
----------

Ultra-Coarse-Graining (UCG) is a more powerful methodology that allows chemical and environmental changes to be captured by modulating the interactions between internal states. In a UCG model, the CG sites can be with a mixed types or "colors" depending on its sorroundings on-the-fly. This algorithm is very important when the same group of atoms of a CG site have different (atomic or electronic) conformations in the simulation system.



The option ``--ucg`` defines the parameters to control the UCG scheme by a comma seperated string with ``key=value`` fields:

  * ``replica=N``, the number of replica to be spawned from a CG frame.
  * ``seed=N``, the random seed for spwaning states.

The option ``--ucg-wf`` defines the weighting functions that's used to calculate the probabilities of states for a CG type. It is also a comma separated string starting with the weighting function name following by a group of ``key=value`` fields:

  * Rapid Local Equilibrium (RLE) Model: ``--ucg-wf RLE,target=MeOH,high=Near,low=Far,rth=4.5,wth=3.5``
  
     * target=[name], the name of the targeted CG type to be spawn.
     * high=[name],low=[name], the names of the CGs types representing sites high or low local densities.
     * rth=[float],wth=[float], the two parameters in RLE model used to calculate local coordination numbers.

.. seealso::
    Jaehyeok Jin and Gregory A. Voth, "Ultra-Coarse-Grained Models Allow for an Accurate and Transferable Treatment of Interfacial Systems", Journal of Chemical Theory and Computation, 14(4), 2180-2197 (2018). doi:10.1021/acs.jctc.7b01173

Examples
--------

*To be added*
