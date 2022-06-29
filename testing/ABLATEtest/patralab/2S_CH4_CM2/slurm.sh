#!/bin/sh
#SBATCH -J ABLATEtest_2S_CH4_CM2
#SBATCH -p mpi #run on mpi queue
#SBATCH --time=0-1:00:00 #day-hour:minute:second
#SBATCH -N 1 #request 1 cores
#SBATCH --mem=2000 #request 2000MG memory
#SBATCH --output=%j.out #your output file
#SBATCH --error=%j.err #your error file
#SBATCH	--mail-type=ALL
#SBATCH	--mail-user=ajbecerr@buffalo.edu
bash 2S_CH4_CM2.sh