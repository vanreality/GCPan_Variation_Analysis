#!/bin/bash
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script=$root/script/callSNP
result=$root/result/mutect2

#samples=$(cat 185samples)
samples=$(head -n 1 185samples)
#samples=$(head -n 11 185samples | tail -n 1)
cd $result || exit
for sample in $samples;
do
  mkdir $sample
  cd $sample || exit
  cp $script/mutect2_filter.slurm ./
  sbatch mutect2.slurm $sample
  cd ..
done
