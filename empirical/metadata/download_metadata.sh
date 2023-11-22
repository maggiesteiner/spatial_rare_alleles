#!/bin/bash

# run from seq_ids directory

# get IDs in WES data
dx run app-swiss-army-knife -y --wait \
    -iin=Rare\ Variants\ in\ Space:/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants,\ PLINK\ format\ -\ final\ release/ukb23158_c8_b0_v1.fam \
    -icmd="awk '{print \$2}' ukb23158_c8_b0_v1.fam > wes_ids.txt"
dx download wes_ids.txt

# get IDs in WGS data
dx run app-swiss-army-knife -y --wait \
    -iin=Rare\ Variants\ in\ Space:/Bulk/Whole\ genome\ sequences/Population\ level\ WGS\ variants,\ PLINK\ format\ -\ interim\ 200k\ release/ukb24305_c8_b0_v1.fam \
    -icmd="awk '{print \$2}' ukb24305_c8_b0_v1.fam > wgs_ids.txt"
dx download wgs_ids.txt

# get metadata

# fields I want
#### participant.eid (ID)
#### participant.p54_i0 (Assessment centre)
#### participant.p21000_i0 (Ethnic background)
#### participant.p1647_i0 (Country of birth UK/elsewhere)
#### participant.p20115_i0 (Country of birth non-UK)
#### participant.p130_i0 (Place of birth in UK - east co-ordinate)
#### participant.p129_i0 (Place of birth in UK - north co-ordinate)

# set dataset environmental variable
dataset="project-GV1P9b8J4p8yz3vQZy9JbJXK:record-GV1v3YjJYq5qqBY4F7Kfgg1J"
# set fields environmental variable
fields="participant.eid,participant.p54_i0,participant.p21000_i0,participant.p1647_i0,participant.p20115_i0,participant.p130_i0,participant.p129_i0"
# generate metadata file
dx extract_dataset "${dataset}" --fields "${fields}" -o participant_metadata.csv


