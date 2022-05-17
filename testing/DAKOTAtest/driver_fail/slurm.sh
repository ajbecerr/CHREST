#!/bin/bash
#SBATCH --partition=debug
#SBATCH --qos=debug
#SBATCH --job-name='chem_sens'
#SBATCH --output=out_chem_sens-%j.out
#SBATCH --error=error_chem_sens-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
module use /projects/academic/chrest/modules
module purge; module load chrest/release
module load dakota/6.15
module load mkl/2020.2
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export DAK_INSTALL=/util/academic/dakota/dakota-6.15.0-release-public-rhel7.x86_64-gui_cli
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH
source $MKL/bin/mklvars.sh
srun dakota -i sensitivity.in

