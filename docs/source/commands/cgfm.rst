cgfm
====

.. include:: common.rst

.. automodule:: mscg.cli.cgfm


* |cli-opt-top-name|

* |cli-opt-traj|

* |cli-opt-model|

* |cli-opt-forcefield|

* For more details about the options ``--ucg`` and ``--ucg-wf``, please read the subsection below for `Ultra Coarse-Graining`_.

Notes
-----

* The script ``cgfm`` doesn't dump the force or coefficient table in the plain-text format. Instead, there are two binary output files (`they are pickle files <https://docs.python.org/3/library/pickle.html>`_),  named by the option ``--save``. By default, they are named *covariance_matrix.p* and *coeffs_matrix.p*.

  * **covariance_matrix.p** stores the covariance matrix for the regression problem before solving it. It can be used as intermediate results for debugging or other purposes.

  * **coeffs_matrix.p** stores the coefficients of B-Spline for each force table after solving the regression problem. Users can use the command `cgdump <cgdump.html>`_ to dump the readable force tables from it.
  
* If this command is `called with in another Python script <../commands.html#call-command-in-python>`_, users may want to access these results directly. In this case, users can pass the value *"return"* to the ``--save`` argument and the coefficient matrix will be returned to the wrapper script as lists. For example::
    
    coeffs = cgfm.main(
        top     = datafile("unary_lj_fluid.top"),
        traj    = datafile("unary_lj_fluid.lammpstrj,frames=20"),
        cut     = 2.50,
        pair    = ['1,1,min=0.9,resolution=0.1'],
        save    = 'return'
    )

* If you have used the old ``MSCGFM`` software, this command is a replacement to the ``newfm`` command in ``MSCGFM``.


Ultra Coarse-Graining
---------------------

Ultra Coarse-Graining (UCG) is a methodology that allows chemical and environmental changes to be captured by modulating the interactions between coarse-grained sites (often as a function of their environment). In a UCG model, the CG sites can comprised of multiple types depending on its surroundings. This algorithm is useful when the group of atoms in a CG site acts differently in different chemical environments. The current implementation primarily support rapid local equilibrium UCG.

The option ``--ucg`` defines the parameters to control the UCG scheme using a comma separated string with the following ``key=value`` fields:

  * ``replica=N``, the number of replica to be spawned from a CG frame.
  * ``seed=N``, the random seed for spawning states.

The option ``--ucg-wf`` defines the weighting functions that's used to calculate the probabilities of states for a CG type. It is also a comma separated string starting with the weighting function name following by a group of ``key=value`` fields. For example:


  * A sample rapid methanol model: ``--ucg-wf RLE,target=MeOH,high=Near,low=Far,rth=4.5,wth=3.5``

     * target=[name], the name of the targeted CG type to spawn.
     * high=[name],low=[name], the names of the CG types representing sites of high or low local densities.
     * rth=[float],wth=[float], the two parameters in the UCG model used to calculate local coordination numbers.

.. seealso::
    Jaehyeok Jin and Gregory A. Voth, "Ultra-Coarse-Grained Models Allow for an Accurate and Transferable Treatment of Interfacial Systems", Journal of Chemical Theory and Computation, 14(4), 2180-2197 (2018). doi:10.1021/acs.jctc.7b01173

Examples
--------

::
    
    cgfm --top ../tests/data/methanol_1728.data \
         --traj ../tests/data/methanol_1728_cg.trr \
         --cut 8.0 --names MeOH --lasso 1.0 \
         --pair model=BSpline,type=MeOH:MeOH,min=2.8,resolution=0.2

::
    
    cgfm --top tests/data/methanol_1728_2s.data \
         --traj tests/data/methanol_1728_2s.trr \
         --cut 8.0 --names CH3,OH \
         --pair model=BSpline type=CH3,CH3 min=2.9 resolution=0.1 order=6 \
         --pair model=BSpline type=CH3,OH min=2.8 resolution=0.1 order=6 \
         --pair model=BSpline type=OH,OH min=2.5 resolution=0.1 order=6 \
         --bond model=BSpline type=CH3,OH min=1.35 max=1.65 resolution=0.01 order=4

::
    
    cgfm --top tests/data/methanol_1728.data \
         --traj run/methanol_slab_cg.trr \
         --cut 8.0 --names MeOH \
         --ucg RLE,type=MeOH
