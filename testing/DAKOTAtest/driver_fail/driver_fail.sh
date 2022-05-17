#!/bin/sh
module use /projects/academic/chrest/modules
module purge; module load chrest/release
module load mkl/2020.2
module load dakota/6.15
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so
export DAK_INSTALL=/util/academic/dakota/dakota-6.15.0-release-public-rhel7.x86_64-gui_cli
export PATH=$DAK_INSTALL/bin:$DAK_INSTALL/share/dakota/test:$PATH
export PYTHONPATH=$DAK_INSTALL/share/dakota/Python:$PYTHONPATH
source $MKL/bin/mklvars.sh
srun dakota -i sensitivity.in