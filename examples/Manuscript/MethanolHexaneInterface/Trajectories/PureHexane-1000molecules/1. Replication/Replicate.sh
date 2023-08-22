#!/bin/bash
#SBATCH --job-name=Rel
#SBATCH --output=Rel-%j.out
#SBATCH --error=Rel-%j.err
#SBATCH --partition=skx-dev
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH --exclusive

module load  intel impi

ibrun /home1/02487/srmani/lammps-7Aug19-RLEUCG/src/lmp_mpi < in.Replicate
