#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --job-name=python

# load environment, e.g. set virtualenv, environment variables, etc


source /scratch/dmoi/miniconda/etc/profile.d/conda.sh
conda activate ML

# Run Jupyter

python coevsankoff_distributed.py
