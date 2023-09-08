#!/bin/bash

lsf_dir=/rsrch5/scratch/biostatistics/wzhang24/CARE/simulation_final_2/pbs
for lsf_file in "$lsf_dir"/sim4_*.lsf; do
	bsub < "$lsf_file"
	sleep 2
done
