#!/bin/bash
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script=$root/script/callSNP
result=$root/result/ref_mutect2
chrom_list=$(cat $root/data/known/chrom_list)

samples=$(cat 185samples)
#samples=$(head -n 160 185samples)
#samples=$(head -n 11 185samples | tail -n 1)
mkdir $result
cd $result || exit
for sample in $samples;
do
  if [ -d $sample ]; then
    continue
  fi
  mkdir $sample
  cd $sample || exit
  cp $script/ref_mutect2.slurm ./
  for i in $chrom_list; do
    sbatch ref_mutect2.slurm $sample $i
  done
  sed -i 's/-L/-XL/g' ref_mutect2.slurm
  sed -i 's/${2}/other/g' ref_mutect2.slurm
  sbatch ref_mutect2.slurm $sample "all"
  cd ..
done
