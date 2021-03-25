#!/bin/bash

#SBATCH -J delly
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
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"
truncated="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSV/truncated_bams"
origin="/lustre/home/acct-clswcc/clswcc-zqd/project/GATK_bam"
restore="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/gatkBams"
delly="/lustre/home/acct-clswcc/clswcc-fsy/software/delly/delly_v0.8.7_linux_x86_64bit"

if [[ $(cat $truncated) =~ $normal ]]; then
	normal_bam="$restore/$normal/$normal.recal.bam"
else
	normal_bam="$origin/$normal/$normal.recal_reads.bam"
fi

if [[ $(cat $truncated) =~ $tumor ]]; then
	tumor_bam="$restore/$tumor/$tumor.recal.bam"
else
	tumor_bam="$origin/$tumor/$tumor.recal_reads.bam"
fi

$delly call -o $normal.$tumor.bcf -g $reference $tumor_bam $normal_bam

printf "%s\tcontrol" $normal >> samples.tsv
printf "%s\ttumor" $tumor >> samples.tsv

$delly filter -p -f somatic -o $normal.$tumor.somatic.bcf -s samples.tsv $normal.$tumor.bcf

bcftools view $normal.$tumor.bcf > $normal.$tumor.vcf
bcftools view $normal.$tumor.somatic.bcf > $normal.$tumor.somatic.vcf