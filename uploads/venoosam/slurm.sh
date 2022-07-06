#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --time=72:00:00
#SBATCH --account=swihart
#SBATCH --nodes=1
#SBATCH --constraint=MRI|NIH
#SBATCH --mem=187000
#SBATCH --job-name="pentane_oxidation"
#SBATCH --output=pentane_oxidation.out
#SBATCH --mail-user=venoosam@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue
#Specifies that the job will be requeued after a node failure.
#The default is that the job will not be requeued.

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR
echo "working directory = "/projects/academic/swihart/RMG/RMG-Py

conda init bash

module load python/py37-anaconda-2020.02
module load rmg/py37

export PATH=~/anaconda3py37/bin:$PATH
source ~/.bashrc

source activate /projects/academic/swihart/RMG/anaconda3py37/envs/rmg_env



eval "$(/util/common/python/py37/anaconda-2020.02/bin/conda shell.bash hook)"

ulimit -s unlimited

source activate rmg_env

#
python rmg.py pentane_oxidation.py
#
echo "All Done!"