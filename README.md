# OpenMSCG

An open-source python package for systematic coarse-graining (including
MSCG/force-matching) in computational chemistry and biology. 

## Quick Guide

**1. Prerequisites**

Install a conda package manager/environment: `Anaconda` or `Miniconda`.

**2. Install package from Anaconda Cloud**

```
conda install -c vothgroup openmscg
```

**4. Test the installation**

```
cginfo
```

## Documentation

See the documentation at https://software.rcc.uchicago.edu/mscg/docs/.

## Notes of Updates & Changes

**0.3.0** (unreleased)

1. Release the feature of substracting tabulated forces from the reference trajectories.
2. In CGIB, add screen output for the ranges of the target variables.

**0.2.2**

1. Enable processing of LAMMPS trajectory in that only part of the system atoms are dumped.
2. Fixed bugs in modeling of multiple types of UCG sites.
3. Allow users to setup number of threads for multithreading acceleration.

**0.2.1**

1. The unit of input values for CGFM angular models is changed from radius to degrees.
2. Add calculation of dihedral angles in CGIB.
3. Add new GaussCut model for the REM method.
4. Change the matrix inversion in CGFM to allow for the singular matrix.
5. Fix a bug in PairList to reduce the need of memory in verlet-list.

**0.1.0**

1. All basic features are released
