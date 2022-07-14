#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --job-name='chem_sens'
#SBATCH --output=out_chem_sens_parallel-%j.out
#SBATCH --error=error_chem_sens_parallel-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
module purge
module load python/py38-anaconda-2020.11
python driver_parallel.py