#!/bin/bash

input_file="files_to_read.txt"

while IFS= read -r line; do
    modified_line=$(echo "$line" | sed 's/\.[^.]*$//;s/\.$//')
    python3 output_moments.py -f "$modified_line"
done < "$input_file"
