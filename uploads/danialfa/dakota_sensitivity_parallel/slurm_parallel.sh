#!/bin/bash
#SBATCH --clusters=faculty
#SBATCH --partition=phoenix
#SBATCH --qos=phoenix
#SBATCH --job-name='test'
#SBATCH --output=out_test-%j.out
#SBATCH --error=error_test-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=100:00:00

#######SBATCH -w, --nodelist=cpn-v11-36

#SBATCH --mail-user=danialfa@buffalo.edu
#SBATCH --mail-type=ALL



module load hdf5/1.12.0-mpi 
module load gcc/10.2.0
module load cmake/3.17.1

export PETSC_DIR=/projects/academic/danialfa/software/petsc
export PETSC_ARCH=arch-linux-c-opt



module load intel/19.5 && \
module load intel-mpi/2019.5 && \
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so && \
module load mkl/2019.5 && \
source $MKL/bin/mklvars.sh intel64 && \
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/academic/danialfa/software/dakota/dependencies/gsl/lib && \
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/academic/danialfa/software/dakota/dependencies/boost/1.49/lib && \
export DAK_INSTALL=/projects/academic/danialfa/software/dakota && \
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH && \
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH


## make it parallel on CCR
export DAKOTA_RUN_PARALLEL=True


srun dakota -i sensitivity_parallel.in

