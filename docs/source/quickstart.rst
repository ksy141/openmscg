Quick Start Guide
=================

.. contents:: :local:

Prerequisites
*************

1. `Python>=3.6 <https://www.python.org/downloads/release/python-360/>`_
2. `PIP>=20.0 <https://pip.pypa.io/en/stable/installing/>`_
3. `numpy>=1.15 <https://pypi.org/project/numpy/>`_
4. `matplotlib>=3.0 <https://pypi.org/project/matplotlib/>`_

Download Precompiled Distributions
**********************************

* For Python 3.6.x users: `wheel-py36 <https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp36-cp36m-linux_x86_64.whl>`_
* For Python 3.7.x users: `wheel-py37 <https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp37-cp37m-linux_x86_64.whl>`_
* For Python 3.8.x users: `wheel-py38 <https://software.rcc.uchicago.edu/download/OpenMSCG-0.1.0-cp38-cp38-linux_x86_64.whl>`_

Quick Install
*************

This package is distributed for all 64-bit Linux platforms. If you're a Microsoft Windows user, please install the simulated Linux environment `cygwin <https://www.cygwin.com/>`_ ::

    pip install OpenMSCG-*.whl


If you don't have the write permission to the folder for Python site-packages, you can install it in your local home directory by using the option ``--user`` ::

    pip install OpenMSCG-*.whl --user

If the package is installed in user's home directory, please export the path for executable scripts to the Linux environment ::

    export PATH=${HOME}/.local/bin:${PATH}


Test the Installtion
********************

Simply run the following command to check the package information ::

    cginfo

Example output ::

                         OpenMSCG Python Package
                        =========================


    > Package Information
    ---------------------

            Name: OpenMSCG
         Version: 0.1.0
         Summary: An open-source package for multi-scale coarse-graining (MSCG) modeling in computational chemistry and
                  biology
       Home-page: https://software.uchicago.edu/mscg/
          Author: The Voth Research Group, The University of Chicago
    Author-email: yuxing@uchicago.edu
         License: GPL
     Description: This software integrates the cutting-edge CG and UCG model generation/simulation algorithms from The
                  Voth Group and provide a user-driven code/data repository for public dissemination of CG/UCG models and
                  parameters. The development is funded by NSF/SSE (No.5-27425), and focused on the broad applicability
                  to the fields of chemistry, biology and the materials sciences in accordance with the National
                  Strategic Supercomputing Initiative (e.g. Principle 1, Objective 2), and the Vision and Strategy for
                  Software for Science, Engineering, and Education (NSF 12-113, Strategies 1-4).
        Platform: linux/x86_64
          Python: >=3.6
        Location: /home/yuxing/.local/lib/python3.6/site-packages/OpenMSCG-0.1.0-py3.6-linux-x86_64.egg/mscg


    > Classes in Package
    --------------------

                BSpline -> mscg.bspline.BSpline
               BondList -> mscg.bondlist.BondList
              CLIParser -> mscg.cli_parser.CLIParser
                 Matrix -> mscg.matrix.Matrix
               PairList -> mscg.pairlist.PairList
                  TIMER -> mscg.timer.TIMER
      TableAngleBSpline -> mscg.tables.angle_bspline.TableAngleBSpline
       TableBondBSpline -> mscg.tables.bond_bspline.TableBondBSpline
       TablePairBSpline -> mscg.tables.pair_bspline.TablePairBSpline
               Topology -> mscg.topology.Topology
             Trajectory -> mscg.trajectory.Trajectory
                 screen -> mscg.verbose.screen
                 tables -> mscg.tables.tables


    Congratulations! The package is successfully installed.
