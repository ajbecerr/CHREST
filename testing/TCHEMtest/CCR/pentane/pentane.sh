#!/bin/sh
export TCHEM_INSTALL_PATH=/projects/academic/chrest/lib/tchem/v2.0.0_25-04-2022_6ae59e8
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