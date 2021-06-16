#!/bin/bash

#SBATCH -J somatic_GWAS
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/snp
res_dir=$root/result/snpGWAS

