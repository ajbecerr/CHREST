
export TCHEM_INSTALL_PATH=${HOME}/Documents/CODE/GitHubCode/install/tchemSSH
exec=$TCHEM_INSTALL_PATH/example/TChem_IgnitionZeroD.x
inputs='../inputs'
this="$exec --chemfile=$inputs/LL2KGB_AllRange_UnitsTesting.yaml \
            --use-cvode=false \
            --samplefile=$inputs/sample.dat \
            --outputfile=IgnSolution.dat \
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
            --ignition-delay-time-file=IgnitionDelayTime.dat \
            --ignition-delay-time-w-threshold-temperature-file=IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"

echo $this
eval $this
