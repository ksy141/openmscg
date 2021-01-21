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