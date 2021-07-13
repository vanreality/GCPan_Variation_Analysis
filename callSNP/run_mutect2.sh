#!/bin/bash
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script=$root/script/callSNP
result=$root/result/mutect2
chrom_list=$(cat $root/data/known/chrom_list)

#samples=$(cat 185samples)
samples=$(head -n 1 185samples)
#samples=$(head -n 11 185samples | tail -n 1)
mkdir $result
cd $result || exit
for sample in $samples;
do
#  if [ -d $sample ]; then
#    continue
#  fi
  mkdir $sample
  cd $sample || exit
  cp $script/mutect2.slurm ./
  for i in $chrom_list; do
    sbatch mutect2.slurm $sample $i
  done
  sed -i 's/-L/-XL/g' mutect2.slurm
  sed -i 's/${2}/other/g' mutect2.slurm
  sbatch mutect2.slurm $sample "all"
  cd ..
done
