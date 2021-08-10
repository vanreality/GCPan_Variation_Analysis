#!/bin/bash

#SBATCH -J pan_delly
#SBATCH -p small
#SBATCH -N 1
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=4

export OMP_NUM_THREADS=4
module load bcftools

normal=$1
tumor=$2
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
reference=$root/data/known/pan.fa
normal_bam=$root/data/panGatkBams/$normal/$normal.recal.bam
tumor_bam=$root/data/panGatkBams/$tumor/$tumor.recal.bam
delly="/lustre/home/acct-clswcc/clswcc-fsy/software/delly/delly_v0.8.7_linux_x86_64bit"
excl="/lustre/home/acct-clswcc/clswcc-fsy/software/delly/delly/excludeTemplates/human.hg38.excl.tsv"

$delly call -x $excl -o $normal.$tumor.bcf -g $reference $tumor_bam $normal_bam

printf "%s\tcontrol\n" $normal >> samples.tsv
printf "%s\ttumor" $tumor >> samples.tsv

$delly filter -p -f somatic -o $normal.$tumor.somatic.bcf -s samples.tsv $normal.$tumor.bcf

bcftools view $normal.$tumor.bcf > $normal.$tumor.vcf
bcftools view $normal.$tumor.somatic.bcf > $normal.$tumor.somatic.vcf