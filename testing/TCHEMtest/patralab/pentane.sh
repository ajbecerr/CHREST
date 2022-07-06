#!/bin/sh
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/gtest/30-03-2022_af29db7/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/kokkos/v3.5.00_30-03-2022_2834f94/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/openblas/release/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/tines/30-03-2022_7dba309/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/chrest/lib/yaml-cpp/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/gcc/11.2.0/gcc/lib64/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/compiler/lib/intel64_lin/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/lib/release/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/patralab/shared/intel/20.2/compilers_and_libraries_2020.2.254/linux/mpi/intel64/libfabric/lib/
export LD_LIBRARY_PATH
export TCHEM_INSTALL_PATH=/cluster/tufts/patralab/shared/chrest/lib/tchem/v2.0.0_25-04-2022_6ae59e8
exec=$TCHEM_INSTALL_PATH/example/TChem_IgnitionZeroD.x
this="$exec --chemfile=inputs/pentane.yaml \
            --use-cvode=false \
            --samplefile=inputs/pentane.dat \
            --outputfile=outputs/IgnSolution.dat \
            --atol-newton=1e-18 \
            --rtol-newton=1e-8\
            --run-constant-pressure=false \
            --max-newton-iterations=20 \
            --tol-time=1e-6 \
            --useYaml=true \
            --dtmax=1e-3 \
            --dtmin=1e-20 \
            --tend=0.1 \
            --time-iterations-per-interval=10 \
            --jacobian-interval=5 \
            --max-time-iterations=5000 \
            --ignition-delay-time-file=outputs/IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=outputs/IgnitionDelayTimeTthreshold.dat \
            --threshold-temperature=1500"
echo $this
eval $this