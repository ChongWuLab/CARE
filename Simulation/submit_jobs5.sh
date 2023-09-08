#!/bin/bash

lsf_dir=/rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/pbs
files=$(find "$lsf_dir" -type f -name "*[1][1-9].lsf" -o -name "*20.lsf")

count=$(echo "$files" | wc -l)
echo "Number of files: $count"
for file in $files
do
    echo "found file: $file"
    bsub < "$file"
    
done