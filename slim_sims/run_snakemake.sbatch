#!/usr/bin/env bash
#SBATCH -J space_snake
#SBATCH --account=pi-jnovembre
#SBATCH --mail-type ALL
#SBATCH --mail-user=mianniravn@uchicago.edu
#SBATCH --partition=jnovembre
#SBATCH --time=12:00:00
#SBATCH --output=snakemake_logs/space_snake_%A_%a.out
#SBATCH --error=snakemake_logs/space_snake_%A_%a.err

module load python
#module load Anaconda3/2022.10
module load R
module load SLiM/3.0
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate snakemake
conda activate ./spatialrarealleles

snakemake --unlock
snakemake --profile slurm --keep-going --rerun-incomplete --printshellcmds
