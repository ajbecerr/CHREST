#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --job-name=ABLATEtest_2S_CH4_CM2
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
module use /projects/academic/chrest/modules
module load chrest/release
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
srun 2S_CH4_CM2.sh