#!/bin/bash

#SBATCH -J extract_vcf
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
res_dir=$root/data/snp
gatk="/lustre/home/acct-clswcc/clswcc-fsy/software/gatk-4.1.9.0/gatk"
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"

# Extract snp from raw vcf files
$gatk SelectVariants -R $reference -V $res_dir/germline_raw.vcf -O $res_dir/germline_snp.vcf --select-type-to-include SNP
$gatk SelectVariants -R $reference -V $res_dir/hanSouth1000genomes/hs.vcf -O $res_dir/hanSouth1000genomes/hs_snp.vcf --select-type-to-include SNP
$gatk SelectVariants -R $reference -V $res_dir/hanBeijing1000genomes/hb.vcf -O $res_dir/hanBeijing1000genomes/hb_snp.vcf --select-type-to-include SNP

# Extract indel from raw vcf files
$gatk SelectVariants -R $reference -V $res_dir/germline_raw.vcf -O $res_dir/germline_indel.vcf --select-type-to-include INDEL
$gatk SelectVariants -R $reference -V $res_dir/hanSouth1000genomes/hs.vcf -O $res_dir/hanSouth1000genomes/hs_indel.vcf --select-type-to-include INDEL
$gatk SelectVariants -R $reference -V $res_dir/hanBeijing1000genomes/hb.vcf -O $res_dir/hanBeijing1000genomes/hb_indel.vcf --select-type-to-include INDEL