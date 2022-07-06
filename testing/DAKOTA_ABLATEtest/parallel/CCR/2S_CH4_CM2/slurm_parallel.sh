#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --job-name='chem_sens'
#SBATCH --output=out_chem_sens_parallel-%j.out
#SBATCH --error=error_chem_sens_parallel-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
module purge
module load dakota/6.15
export DAK_INSTALL=/util/academic/dakota/dakota-6.15.0-release-public-rhel7.x86_64-gui_cli
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH
export DAKOTA_RUN_PARALLEL=True
srun dakota -i sensitivity_parallel.in