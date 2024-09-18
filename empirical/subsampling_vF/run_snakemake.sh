#!/usr/bin/env bash
#SBATCH -J ukb_sfs
#SBATCH --account=pi-jnovembre
#SBATCH --mail-type ALL
#SBATCH --mail-user=steinerm@uchicago.edu
#SBATCH --partition=broadwl
#SBATCH --time=02:00:00

module load python
module load R
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate snakemake

snakemake --profile slurm
   

