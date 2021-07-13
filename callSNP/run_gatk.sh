#!/bin/bash

script="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSNP"
result="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/panGatkBams"

samples=$(cat 185samples)
#samples=$(head -n 1 185samples)
#samples=$(head -n 11 185samples | tail -n 1)
mkdir $result
cd $result || exit
for sample in $samples;
do
  normal=$(echo $sample | awk -F "." '{print $1}')
  tumor=$(echo $sample | awk -F "." '{print $2}')
  if [ -d $normal ] && [ -d $tumor ]; then
    continue
  fi
  mkdir $normal
  cd $normal || exit
  cp $script/gatk.slurm ./
  sbatch gatk.slurm $normal
  cd ..

  mkdir $tumor
  cd $tumor || exit
  cp $script/gatk.slurm ./
  sbatch gatk.slurm $tumor
  cd ..
done
