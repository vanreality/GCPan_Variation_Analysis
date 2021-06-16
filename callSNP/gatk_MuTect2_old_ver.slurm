#!/bin/bash

#SBATCH -J merge
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

gatk="/lustre/home/acct-clswcc/clswcc-fsy/software/gatk-4.1.9.0/gatk"
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
samples=$(cat $root/script/callSV/185samples)
vcf_dir="/lustre/home/acct-clswcc/clswcc-zqd/project/VariantCalling/MuTect2/vcf"
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"
knownsites_dir="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/known"

for sample in $samples ; do
    gzip -d $vcf_dir/$sample.noContEst.vcf.gz

done