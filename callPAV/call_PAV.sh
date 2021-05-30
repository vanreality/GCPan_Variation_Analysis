#!/bin/bash

#samples=$(head -n 1 185samples)
samples=$(cat 185samples)
for sample in $samples;
do

#batch process

normal=$(echo $sample | awk -F "." '{print $1}')
tumor=$(echo $sample | awk -F "." '{print $2}')

#sbatch cov_ref.slurm $normal "Normal"
#sbatch cov_ref.slurm $tumor "Tumor"

sbatch cov_pan.slurm $normal "Normal"
sbatch cov_pan.slurm $tumor "Tumor"

done;
