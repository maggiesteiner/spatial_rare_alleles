#!/bin/bash

#dx run app-swiss-army-knife -y --wait \
#    -iin=Rare\ Variants\ in\ Space:/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ final\ release/ukb23158_c8_b0_v1.fam \
#    -icmd="awk '{print \$2}' ukb23158_c8_b0_v1.fam > wes_ids_all.txt"
#    dx mv wes_ids_all.txt maggie_pipeline_v3
#    dx download maggie_pipeline_v3/wes_ids_all.txt
#    mv wes_ids_all.txt data/metadata

dataset="project-GV1P9b8J4p8yz3vQZy9JbJXK:record-GV1v3YjJYq5qqBY4F7Kfgg1J"

### FIELD DESCRIPTORS ###
### participant.eid: participant ID
### participant.p22001: genetic sex
### participant.p21022: age at recruitment
### participant.p12144_i2: height

fields="participant.eid,participant.p21022,participant.p22001,participant.p12144_i2"
dx extract_dataset "${dataset}" --fields "${fields}" -o participant_metadata_gwas.csv
mv participant_metadata_gwas.csv data/metadata
