#!/bin/bash

#SBATCH -J pan_svaba
#SBATCH -p small
#SBATCH -N 1
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=16

normal=$1
tumor=$2
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
known_indel="/lustre/home/acct-clswcc/clswcc-zqd/data/known/dbsnp_indel.vcf"
reference=$root/data/known/pan.fa
normal_bam=$root/data/panGatkBams/$normal/$normal.recal.bam
tumor_bam=$root/data/panGatkBams/$tumor/$tumor.recal.bam
svaba="/lustre/home/acct-clswcc/clswcc-fsy/software/svaba/bin/svaba"

$svaba run -t $tumor_bam -n $normal_bam -p 16 -D $known_indel -a $normal.$tumor -G $reference