#!/bin/sh
#SBATCH -J DAKOTA_TCHEMtest_pentane
#SBATCH -p batch #run on batch queue
#SBATCH --time=0-1:00:00 #day-hour:minute:second
#SBATCH -N 1 #request 1 cores
#SBATCH --mem=2000 #request 2000MG memory
#SBATCH --output=outfile #your output file
#SBATCH --error=errfile #your error file
module load dakota/6.13
export DAK_INSTALL=/cluster/tufts/hpc/tools/dakota/6.13.1/dakota
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH
srun dakota -i sensitivity.in