#!/bin/sh
module purge
module load lapack/3.5.0
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/ablate/v0.8.14_19-05-2022_acf63a8/src/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/petsc/v3.17.1_19-05-2022_673874f/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/gcc/11.2.0/gcc/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/lib/release/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/libfabric/lib/
export LD_LIBRARY_PATH
/cluster/tufts/patralab/shared/chrest/lib/ablate/v0.8.14_19-05-2022_acf63a8/ablate --input ignitionDelay2S_CH4_CM2.yaml