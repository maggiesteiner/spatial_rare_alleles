#!/usr/bin/env bash
#SBATCH -J ukb_pca
#SBATCH --mem=50G
#SBATCH --time=4:00:00

ml python
source /software/python-anaconda-2020.02-el7-x86_64/etc/profile.d/conda.sh
conda activate ukb_pca
ml plink

#### STEP 1: Merge files across autosomes ####

cd plink_files
mkdir -p merged_files
for chrom in {1..22}; do
    echo "Chromosome ${chrom}"
    base_file="c${chrom}_filtered"
    if [ "$chrom" -eq 1 ]; then
        plink --file "$base_file" --make-bed --out "merged_files/merged_dataset"
    else
        plink --bfile "merged_files/merged_dataset" --bmerge "$base_file" --allow-no-sex --make-bed --out "merged_files/merged_dataset_temp"
        plink --bfile "merged_files/merged_dataset_temp" --make-bed --out "merged_files/merged_dataset"
        rm -f "merged_files/merged_dataset_temp"
    fi
done






