#!/bin/bash

# note - log in to dnanexus before running! and go to smartpca_genotypes directory

path=$"Rare\ Variants\ in\ Space:/smartpca_genotypes/"

dx run app-swiss-army-knife -y \
    --instance-type '{"main": "mem1_ssd1_v2_x36","*": "mem1_ssd1_v2_x36"}' \
    -iin=merged_filtered.bed \
    -iin=merged_filtered.bim \
    -iin=merged_filtered.fam \
    -iin=used_in_pca.tsv \
    -icmd="plink2 --bfile merged_filtered --keep used_in_pca.tsv --pca approx 20 --out ukb_pca_plink2"
