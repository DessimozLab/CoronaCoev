#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --job-name=python
#SBATCH --output=master_slurm.out


source /work/FAC/FBM/DBC/cdessim2/default/dmoi/condaenvs/etc/profile.d/conda.sh

conda activate ML2
# Run Jupyter
python coevsankoff_distributed_delayed_remoteopen.py
