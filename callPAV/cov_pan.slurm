#!/bin/bash

#SBATCH -J cov_pan
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=cov_pan.out
#SBATCH --error=cov_pan.err

sample=$1
sampleType=$2
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
result=$root/data/refPAV/pan
gene_annotation=$root/data/refPAV/pan.pTpG.gff3
bam=/lustre/home/acct-clswcc/clswcc-zqd/project/Pangenome/Map2Pan/$sampleType/$sample/$sample.bam

bedtools genomecov -bga -split -ibam $bam > $result/$sample.bed

python $root/script/callPAV/gene_coverage_v0.py $result/$sample.bed $gene_annotation 0 0.8 $sample > $result/cov/$sample.cov

cat $result/cov/$sample.cov | grep -v "#" >> $result/$sampleType.pan.cov

rm $result/$sample.bed
