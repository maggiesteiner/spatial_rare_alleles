#!/usr/bin/env bash
#SBATCH -J convertf
#SBATCH --mem=50G
#SBATCH --time=4:00:00

module load python
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate ukb_pca

awk '{for(i=6;i<=NF;i++) if($i == -9) $i=1}1' plink_files/merged_filtered.fam > tmp.fam && mv tmp.fam plink_files/merged_filtered.fam
convertf -p paramfile.convertf
