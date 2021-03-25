#!/bin/bash

#SBATCH -J svaba
#SBATCH -p small
#SBATCH -N 1
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=16

normal=$1
tumor=$2
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"
known_indel="/lustre/home/acct-clswcc/clswcc-zqd/data/known/dbsnp_indel.vcf"
truncated="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSV/truncated_bams"
origin="/lustre/home/acct-clswcc/clswcc-zqd/project/GATK_bam"
restore="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/gatkBams"
svaba="/lustre/home/acct-clswcc/clswcc-fsy/software/svaba/bin/svaba"

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

$svaba run -t $tumor_bam -n $normal_bam -p 16 -D $known_indel -a $normal.$tumor -G $reference

