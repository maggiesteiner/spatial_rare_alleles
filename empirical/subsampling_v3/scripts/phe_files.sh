#!/bin/bash

# Input CSV file
input_csv="data/metadata/participant_metadata_gwas.csv"

# Output files
output_tsv1="GWAS_input/covar.phe"
output_tsv2="GWAS_input/height.phe"

# Create the first TSV file with desired columns
echo -e "FID\tIID\tage_at_recruitment\tgenetic_sex" > "$output_tsv1"

# Process CSV file for the first output and append the data
awk -F',' 'BEGIN {OFS="\t"} NR>1 {for(i=1; i<=NF; i++) if ($i == "") $i = "NA"; print $1, $1, $2, $3}' "$input_csv" >> "$output_tsv1"

# Create the second TSV file with desired columns
echo -e "FID\tIID\theight" > "$output_tsv2"

# Process CSV file for the second output and append the data
awk -F',' 'BEGIN {OFS="\t"} NR>1 {for(i=1; i<=NF; i++) if ($i == "") $i = "NA"; print $1, $1, $4}' "$input_csv" >> "$output_tsv2"
