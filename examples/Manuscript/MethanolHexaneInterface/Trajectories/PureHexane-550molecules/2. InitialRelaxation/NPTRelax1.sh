#!/bin/bash
#SBATCH --job-name=Rel
#SBATCH --output=Rel-%j.out
#SBATCH --error=Rel-%j.err
#SBATCH --partition=skx-dev
#SBATCH --time=02:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH --exclusive

module load  intel impi

ibrun /home1/02487/srmani/lammps-7Aug19-RLEUCG/src/lmp_mpi < in.NPTRelax1 > lammps1.out
