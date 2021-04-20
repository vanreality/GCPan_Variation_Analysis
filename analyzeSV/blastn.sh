#!/bin/bash

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzeSV
manta_dir=$root/result/manta
result_dir=$root/result/blastn
if [ ! -d $result_dir ]; then
    mkdir $result_dir
fi

#samples=$(head $script_dir/185samples -n 1)
#samples=$(head $script_dir/185samples -n 10)
samples=$(cat $script_dir/185samples)
for sample in $samples; do
  if [ -d $result_dir/$sample ]; then
      continue
  fi

  mkdir $result_dir/$sample
  vcf_dir=$manta_dir/$sample/result/results/variants
  if [ -f $vcf_dir/germlineSV.vcf.fasta ]; then
    cd $result_dir/$sample || exit
    cp $root/script/analyzeSV/blastn.slurm ./
    sbatch blastn.slurm $vcf_dir/germlineSV.vcf.fasta ./germlineSV.novelSeq.blastn.out
  fi

  if [ -f $vcf_dir/somaticSV.vcf.fasta ]; then
    cd $result_dir/$sample || exit
    cp $root/script/analyzeSV/blastn.slurm ./
    sbatch blastn.slurm $vcf_dir/somaticSV.vcf.fasta ./somaticSV.novelSeq.blastn.out
  fi
done;