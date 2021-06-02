#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=400G
#SBATCH --time=72:00:00
#SBATCH --job-name=fast-tree

/scratch/dmoi/software/veryfasttree/VeryFastTree -threads 30 -fastexp 3 -ext AVX2 -gtr -nt /scratch/dmoi/datasets/covid_data/msa_0501/msa_0501.fasta.cleaned.fasta  > /scratch/dmoi/datasets/covid_data/msa_0501/msa_0501_tree.nwk
