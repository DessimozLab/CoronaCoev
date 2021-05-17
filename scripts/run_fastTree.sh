#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --job-name=fast-tree

/scratch/dmoi/datasets/covid_data/msa_0501/FastTreeDbl -nt -gtr < /scratch/dmoi/datasets/covid_data/msa_0501/msa_0501.fasta.cleaned.fasta > /scratch/dmoi/datasets/covid_data/msa_0501/bigtree.newick
