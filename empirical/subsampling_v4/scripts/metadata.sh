#!/bin/bash

dx run app-swiss-army-knife -y --wait \
    -iin=Rare\ Variants\ in\ Space:/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ final\ release/ukb23158_c8_b0_v1.fam \
    -icmd="awk '{print \$2}' ukb23158_c8_b0_v1.fam > wes_ids_all.txt"
    dx mv wes_ids_all.txt maggie_pipeline_v3
    dx download maggie_pipeline_v3/wes_ids_all.txt
    mv wes_ids_all.txt data/metadata

dataset="project-GV1P9b8J4p8yz3vQZy9JbJXK:record-GV1v3YjJYq5qqBY4F7Kfgg1J"

### FIELD DESCRIPTORS ###
### participant.eid: participant ID
### participant.p22006: genetic ethnic grouping (White British)
### participant.p22020: used in PCA (QC + relatedness filtering)
### participant.p54_i0: assessment centre
### participant.p130_i0: place of birth East coordinate (in UK)
### participant.p129_i0: place of birth North coordinate (in UK)
### participant.p1647_i0: country of birth (UK/elsewhere)

fields="participant.eid,participant.p22006,participant.p22020,participant.p54_i0,participant.p130_i0,participant.p129_i0,participant.p1647_i0"
dx extract_dataset "${dataset}" --fields "${fields}" -o participant_metadata_all.csv
mv participant_metadata_all.csv data/metadata
