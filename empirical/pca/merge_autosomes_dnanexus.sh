#!/bin/bash

# note - log in to dnanexus before running! and go to smartpca_genotypes directory

path=$"Rare\ Variants\ in\ Space:/smartpca_genotypes/"

dx run app-swiss-army-knife -y\
    -iin="${path}c1_filtered.bed" \
    -iin="${path}c1_filtered.bim" \
    -iin="${path}c1_filtered.fam" \
    -iin="${path}c2_filtered.bed"    \
    -iin="${path}c2_filtered.bim"    \
    -iin="${path}c2_filtered.fam"    \
    -iin="${path}c3_filtered.bed"    \
    -iin="${path}c3_filtered.bim"    \
    -iin="${path}c3_filtered.fam"    \
    -iin="${path}c4_filtered.bed"    \
    -iin="${path}c4_filtered.bim"    \
    -iin="${path}c4_filtered.fam"    \
    -iin="${path}c5_filtered.bed"    \
    -iin="${path}c5_filtered.bim"    \
    -iin="${path}c5_filtered.fam"    \
    -iin="${path}c6_filtered.bed"    \
    -iin="${path}c6_filtered.bim"    \
    -iin="${path}c6_filtered.fam"    \
    -iin="${path}c7_filtered.bed"    \
    -iin="${path}c7_filtered.bim"    \
    -iin="${path}c7_filtered.fam"    \
    -iin="${path}c8_filtered.bed"    \
    -iin="${path}c8_filtered.bim"    \
    -iin="${path}c8_filtered.fam"    \
    -iin="${path}c9_filtered.bed"    \
    -iin="${path}c9_filtered.bim"    \
    -iin="${path}c9_filtered.fam"    \
    -iin="${path}c10_filtered.bed"    \
    -iin="${path}c10_filtered.bim"    \
    -iin="${path}c10_filtered.fam"    \
    -iin="${path}c11_filtered.bed"    \
    -iin="${path}c11_filtered.bim"    \
    -iin="${path}c11_filtered.fam"    \
    -iin="${path}c12_filtered.bed"    \
    -iin="${path}c12_filtered.bim"    \
    -iin="${path}c12_filtered.fam"    \
    -iin="${path}c13_filtered.bed"    \
    -iin="${path}c13_filtered.bim"    \
    -iin="${path}c13_filtered.fam"    \
    -iin="${path}c14_filtered.bed"    \
    -iin="${path}c14_filtered.bim"    \
    -iin="${path}c14_filtered.fam"    \
    -iin="${path}c15_filtered.bed"    \
    -iin="${path}c15_filtered.bim"    \
    -iin="${path}c15_filtered.fam"    \
    -iin="${path}c16_filtered.bed"    \
    -iin="${path}c16_filtered.bim"    \
    -iin="${path}c16_filtered.fam"    \
    -iin="${path}c17_filtered.bed"    \
    -iin="${path}c17_filtered.bim"    \
    -iin="${path}c17_filtered.fam"    \
    -iin="${path}c18_filtered.bed"    \
    -iin="${path}c18_filtered.bim"    \
    -iin="${path}c18_filtered.fam"    \
    -iin="${path}c19_filtered.bed"    \
    -iin="${path}c19_filtered.bim"    \
    -iin="${path}c19_filtered.fam"    \
    -iin="${path}c20_filtered.bed"    \
    -iin="${path}c20_filtered.bim"    \
    -iin="${path}c20_filtered.fam"    \
    -iin="${path}c21_filtered.bed"    \
    -iin="${path}c21_filtered.bim"    \
    -iin="${path}c21_filtered.fam"    \
    -iin="${path}c22_filtered.bed"    \
    -iin="${path}c22_filtered.bim"    \
    -iin="${path}c22_filtered.fam"    \
    -iin="${path}used_in_pca.tsv" \
    -iin="${path}all_files.txt" \
    -icmd="plink -bfile c1_filtered --merge-list all_files.txt --make-bed --out merged_filtered"

