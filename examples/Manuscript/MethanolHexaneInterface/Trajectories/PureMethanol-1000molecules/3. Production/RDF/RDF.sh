#!/bin/bash
#SBATCH --job-name=Rel
#SBATCH --output=Rel-%j.out
#SBATCH --error=Rel-%j.err
#SBATCH --partition=development
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH --exclusive

module load  intel impi

Start=0
End=5000000
Dump=5000000
N_Iter=`echo $Start $End $Dump | awk '{print ($2-$1)/$3}'`
echo ${N_Iter}
for i in `seq 1 1 ${N_Iter}`
do
CStart=`echo $Start $i $Dump | awk '{print $1+($2-1)*$3}'`
CEnd=`echo $CStart $Dump | awk '{print $1+$2}'`
NFiles=`echo $Dump 1000 | awk '{print $1/$2}'`
Time=${CStart}_${CEnd}
echo ${CStart} ${CEnd} ${NFiles} ${Time}
cp in.RDF in.RDF_${Time}
sed -i "s/DUMP/$Dump/"g in.RDF_${Time}
sed -i "s/NFILES/$NFiles/"g in.RDF_${Time}
sed -i "s/START/$CStart/"g in.RDF_${Time}
sed -i "s/END/$CEnd/"g in.RDF_${Time}
sed -i "s/TIME/$Time/"g in.RDF_${Time}
ibrun /home1/02487/srmani/lammps-7Aug19-RLEUCG/src/lmp_mpi < in.RDF_${Time}
done
