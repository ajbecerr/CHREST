#!/bin/sh
#SBATCH -J ABLATEtest_2S_CH4_CM2
#SBATCH -p batch #run on batch queue
#SBATCH -N 1 #request 1 cores
#SBATCH --mem=2000 #request 2000MG memory
#SBATCH --output=outfile #your output file
#SBATCH --error=errfile #your error file
#SBATCH	--mail-type=ALL
#SBATCH	--mail-user=ajbecerr@buffalo.edu
srun 2S_CH4_CM2.sh