Quick Start Guide
=================

.. contents:: :local:

Prerequisites
-------------

``Conda`` package and environment management tool, `Anaconda <https://www.anaconda.com/products/individual>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`.

Quick Install
-------------

With `conda` setup in the environment, install the software simply by::

    conda install -c vothgroup openmscg

Test the Installtion
--------------------

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

