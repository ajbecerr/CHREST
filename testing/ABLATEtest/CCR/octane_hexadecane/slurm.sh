#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --job-name=ABLATEtest_octane_hexadecane
#SBATCH --output=%.out
#SBATCH --error=%.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=72:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
module use /projects/academic/chrest/modules
module load chrest/release
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun octane_hexadecane.sh