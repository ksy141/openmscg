Example of 1-Site Water Box by Relative Entropy Method
======================================================

0. Download the reference trajectory
```
   wget --no-check-certificate https://software.rcc.uchicago.edu/mscg/downloads/openmscg-examples/reference.lammpstrj
```

1. Prepare the LAMMPS input for the testing system

  * system.data -- LAMMPS data-file
  * in.lmp      -- Input script for LAMMPS
  * md.inp      -- Shell script to run LAMMPS

2. Generate reference data

  * reference.table     -- Reference potential table (target)
  * reference.lammpstrj -- Trjectory generated from the reference table
  * gen_ref.sh          -- MSCG(cgderiv) script to create reference parameters
  * model_ref.p         -- Generated reference model

3. Test the setup of MSCG/cgderiv for a single iteration

  * cgderiv.sh     -- Script for single iteration
  * dump.lammpstrj -- Trial trajectory for testing

4. Construct initial model by IB

  * cgib.sh   -- calculate initial potential model from Boltzmann-Inversion
  * model.txt -- Initial model

5. Run iterative aapproach

  * cgrem.sh         -- MSCG(cgrem) script
  * Pair_SL-SL.table -- LAMMPS table of model potential after iterations
