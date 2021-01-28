Performance
===========

This section introduce how to run the code parallelly to improve the performance on a multi-core CPU or cluster platforms.


Multithreading
--------------

OpenMSCG utilizes the Intel MKL library, which provides high-performance multithreading for linear algebra operations. Users can tune this feature via the environment variable `OPENMSCG_THREADS`. The value of this variable will decide how many threads/cores to be used by OpenMSCG.

The command below will allow OpenMSCG to use four threads::
    
    $ export OPENMSCG_THREADS=4

The command below will limit OpenMSCG in a single-thread mode::
    
    $ export OPENMSCG_THREADS=1

Setting its value to 0 or not setting it will allow OpenMSCG & MKL determine the number of available cores on the maachine automatically and lanuch the same number of threads. 


Distributed Workflow
--------------------

To speed up the modeling work in a cluster platform, the best practice is to process diferent segments of trajectories on different computing nodes independently. After distributed processing, the results can be combined together to generate the final models. An example of distributed force-matching work is described as below.

The reference trajectory is splitted in 5 files: `traj1.dcd, traj2.dcd, ... traj5.dcd`. The following batch script can be used on a cluster using SLURM the scheduling system.

.. code-block::
    :caption: job.sbatch
    
    #!/bin/sh
    
    #SBATCH --nodes=1
    #SBATCH --exclusive
    #SBATCH --array=1-5
    
    cgfm @cgfm.args --traj traj${SLURM_ARRAY_TASK_ID}.dcd --save result/${SLURM_ARRAY_TASK_ID}
    
In this example, an array of 5 jobs is dispatched, each of which uses all cores of a node to process one trajectory file. The common arguments for all sub-tasks are provided in the file `cgfm.args`, such as `--top` and `--model`, while the argument with different values, i.e., `--traj`, is specified in the command line. The environment variable `$SLURM_ARRAY_TASK_ID` is assigned by SLURM for each job to process different file and save the result in a different checkpoint file in the folder `result`.

After the processing of each file is finished. The following command can be used to read all results and recalculate the model from the summations of them::
    
    cgnes --equation result/*

Besides processing different files, another approach is to process different parts in the same file, with the following tip::
    
    --traj full.dcd,skip=$((SLURM_ARRAY_TASK_ID * 200)),frames=200

Supposing the trajectory `traj.dcd` contains 1000 frames, the argument above allows each sub job in an array to process every 200 frames of it.

