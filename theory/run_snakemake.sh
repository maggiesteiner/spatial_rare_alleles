#!/usr/bin/env bash
#SBATCH -J rvsfs
#SBATCH --mem=10000
#SBATCH --time=4-00:00:00
#SBATCH --account=pi-jnovembre
#SBATCH --mail-type ALL
#SBATCH --mail-user=steinerm@uchicago.edu
#SBATCH --partition=jnovembre

module load python
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate snakemake

snakemake --cores 14 --nolock 

