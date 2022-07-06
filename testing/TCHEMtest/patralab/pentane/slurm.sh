#!/bin/sh
#SBATCH -J TCHEMtest_pentane
#SBATCH -p batch #run on batch queue
#SBATCH --time=0-1:00:00 #day-hour:minute:second
#SBATCH -N 1 #request 1 cores
#SBATCH --mem=2000 #request 2000MG memory
#SBATCH --output=outfile #your output file
#SBATCH --error=errfile #your error file
bash pentane.sh