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
module load python/py37-anaconda-2020.02
module load rmg/py37
eval "$(/util/common/python/py37/anaconda-2020.02/bin/conda shell.bash hook)"
conda activate rmg_env
python /projects/academic/swihart/Venus/computations/RMG-jobs/job1/RMG-Py/rmg.py heptane.py