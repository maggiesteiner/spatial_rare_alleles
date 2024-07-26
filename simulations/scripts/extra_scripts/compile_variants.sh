#!/bin/bash

# Extract the header from the first variant file and create the output file
header_extracted=false
for file in *.variants.tsv; do
    if [ "$header_extracted" = false ]; then
        head -n 1 "$file" > all_variants.tsv
        header_extracted=true
    fi
done

# Concatenate the contents of all other variant files (excluding headers)
for file in *.variants.tsv; do
    tail -n +2 "$file" >> all_variants.tsv
done

echo "Concatenation completed. Output saved as all_variants.tsv"

