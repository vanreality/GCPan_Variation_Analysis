#!/bin/bash

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzeSV
blastn_dir=$root/result/blastn
manta_dir=$root/result/manta

#samples=$(head $script_dir/185samples -n 1)
#samples=$(head $script_dir/185samples -n 10)
samples=$(cat $script_dir/185samples)

for sample in $samples ; do
    cd $blastn_dir/$sample || exit
    cp $script_dir/filter_blastn.slurm ./
    cp $script_dir/filter_blastn.py ./
    if [ -f germlineSV.blastn.out ]; then
      if [ ! -f germlineSV.blastn.out.filtered ]; then
        sbatch blastn_filter.slurm germlineSV.blastn.out $manta_dir/$sample/result/results/variants/germlineSV.vcf.fasta
      fi
    fi
    if [ -f somaticSV.blastn.out ]; then
      if [ ! -f somaticSV.blastn.out.filtered ]; then
        sbatch blastn_filter.slurm somaticSV.blastn.out $manta_dir/$sample/result/results/variants/somaticSV.vcf.fasta
      fi
    fi
done