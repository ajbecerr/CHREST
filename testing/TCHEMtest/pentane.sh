#!/bin/sh
#SBATCH --constraint=CPU-Gold-6130
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
module load chrest/release
export TCHEM_INSTALL_PATH=${HOME}/Documents/CODE/GitHubCode/install/tchemSSH
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
            --ignition-delay-time-w-threshold-temperature-file=outputs/IgnitionDelayTimeTthreshold.dat
            --threshold-temperature=1500"

echo $this
eval $this