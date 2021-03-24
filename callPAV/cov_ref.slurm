#!/bin/bash

#SBATCH -J cov_ref
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=cov_ref.out
#SBATCH --error=cov_ref.err

sample=$1
sampleType=$2
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
origin="/lustre/home/acct-clswcc/clswcc-zqd/project/GATK_bam"
restore=$root/data/gatkBams
result=$root/data/refPAV/ref
gene_annotation=$root/data/refPAV/gencode.v30.annotation.pTpG.gff3
truncated=$root/script/callPAV/truncated_bams

if [[ `cat $truncated` =~ $sample ]]; then
	bam="$restore/$sample/$sample.recal.bam"
else
	bam="$origin/$sample/$sample.recal_reads.bam"
fi

bedtools genomecov -bga -split -ibam $bam > $result/$sample.bed

python $root/script/callPAV/gene_coverage_v0.py $result/$sample.bed $gene_annotation 0 0.8 $sample > $result/cov/$sample.cov

cat $result/cov/$sample.cov | grep -v "#" >> $result/$sampleType.ref.cov

rm $result/$sample.bed
