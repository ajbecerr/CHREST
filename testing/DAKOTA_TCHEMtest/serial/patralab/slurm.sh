#!/bin/sh
#SBATCH -J DAKOTA_TCHEMtest
#SBATCH -N 1 # nodes requested
#SBATCH --mem=2000 # memory in Mb
#SBATCH -o outfile # send stdout to outfile
#SBATCH -e errfile # send stderr to errfile
#SBATCH	--mail-type=ALL
#SBATCH	--mail-user=ajbecerr@buffalo.edu
module purge
module load dakota/6.13
export DAK_INSTALL=/cluster/tufts/hpc/tools/dakota/6.13.1/dakota
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH
srun dakota -i sensitivity.in