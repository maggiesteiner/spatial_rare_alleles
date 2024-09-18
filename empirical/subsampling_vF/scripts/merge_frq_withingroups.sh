#!/bin/bash

# Directory containing the TSV files
input_dir="data/frq"

# Loop over chromosomes 1 to 22
for chrom in {1..22}; do
    # Define the output file name for each chromosome
    output_file="${input_dir}/merged_chr${chrom}_withingroups_all.tsv"
    
    # Define the four specific TSV files for the current chromosome
    tsv_files=(
        "${input_dir}/chr${chrom}_within1epsilon_pca_all.frq.strat"
        "${input_dir}/chr${chrom}_notwithin5epsilon_geo_all.frq.strat"
        "${input_dir}/chr${chrom}_notwithin1epsilon_pca_all.frq.strat"
        "${input_dir}/chr${chrom}_within5epsilon_geo_all.frq.strat"
    )
    
    # Merge the files, skipping the header for all except the first file
    awk 'FNR==1 && NR!=1{next;}{print}' "${tsv_files[@]}" > $output_file

    echo "Merged files for chromosome ${chrom} into ${output_file}"
done

