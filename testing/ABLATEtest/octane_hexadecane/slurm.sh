#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --job-name='octane_hexadecane'
#SBATCH --output=octane_hexadecane-%j.out
#SBATCH --error=octane_hexadecane-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
bash octane_hexadecane.sh