#!/bin/bash

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzeSV
blastn_dir=$root/result/blastn

#samples=$(head $script_dir/185samples -n 1)
#samples=$(head $script_dir/185samples -n 10)
samples=$(cat $script_dir/185samples)

for sample in $samples ; do
    cd $blastn_dir/$sample || exit
    cp $script_dir/filter_blastn.slurm ./
    cp $script_dir/filter_blastn.py ./
    if [ -f germlineSV.blastn.out ]; then
      if [ ! -f germlineSV.blastn.out.filtered ]; then
        sbatch filter_blastn.slurm germlineSV.blastn.out
      fi
    fi
    if [ -f somaticSV.blastn.out ]; then
      if [ ! -f somaticSV.blastn.out.filtered ]; then
        sbatch filter_blastn.slurm somaticSV.blastn.out
      fi
    fi
done