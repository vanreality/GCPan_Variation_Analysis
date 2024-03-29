#!/bin/bash

#SBATCH -J manta
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

normal=$1
tumor=$2
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"
manta="/lustre/home/acct-clswcc/clswcc-fsy/miniconda3/envs/mypy_py27/bin/configManta.py"
truncated="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSV/truncated_bams"
origin="/lustre/home/acct-clswcc/clswcc-zqd/project/GATK_bam"
restore="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/gatkBams"

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

conda activate mypy_py27
$manta --normalBam $normal_bam --tumorBam $tumor_bam --referenceFasta $reference --runDir result
$PWD/result/runWorkflow.py -m local -j 16

cd $PWD/result/results/variants || exit
gzip -d *.vcf.gz
mv diploidSV.vcf germlineSV.vcf






