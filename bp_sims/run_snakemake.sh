#!/usr/bin/env bash
#SBATCH -J bp_sims
#SBATCH --account=pi-jnovembre
#SBATCH --partition=broadwl
#SBATCH --time=08:00:00

module load python
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate snakemake

snakemake --profile slurm
   

