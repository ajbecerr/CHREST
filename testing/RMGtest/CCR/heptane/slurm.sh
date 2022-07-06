#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --time=72:00:00
#SBATCH --account=swihart
#SBATCH --nodes=1
#SBATCH --constraint=MRI|NIH
#SBATCH --mem=187000
#SBATCH --job-name="RMGtest_heptane"
#SBATCH --output=heptane.out
#SBATCH --error=heptane.err
#SBATCH --mail-user=ajbecerr@buffalo.edu
#SBATCH --mail-type=ALL
eval "$(/util/common/python/py37/anaconda-2020.02/bin/conda shell.bash hook)"
conda activate rmg_env
pip install numpy==1.16.1
export PATH=/util/common/rmg/3.0.0/RMG-Py:$PATH
export PYTHONPATH=/util/common/rmg/3.0.0/RMG-Py:$PYTHONPATH
python /util/common/rmg/3.0.0/RMG-Py/rmg.py heptane.py