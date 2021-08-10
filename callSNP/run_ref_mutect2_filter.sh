#!/bin/bash
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script=$root/script/callSNP
result=$root/result/ref_mutect2

#samples=$(cat 185samples)
samples=$(head -n 1 185samples)
#samples=$(head -n 11 185samples | tail -n 1)
cd $result || exit
for sample in $samples;
do
  mkdir $sample
  cd $sample || exit
  cp $script/ref_mutect2_filter.slurm ./
  sbatch ref_mutect2_filter.slurm $sample
  cd ..
done
