#!/bin/bash

#SBATCH -J ref_mutect2
#SBATCH -p arm128c256g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=20

# arm software
module use /lustre/share/spack/modules/kunpeng920/linux-centos7-aarch64
module load gatk/4.2.0.0-gcc-9.3.0-openblas
module load bwa/0.7.17-gcc-9.3.0
module load gcc

# data
sample=$1
normal=$(echo $sample | awk -F "." '{print $1}')
tumor=$(echo $sample | awk -F "." '{print $2}')
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
reference=$root/data/known/hg38.fa
pon=$root/data/known/somatic-hg38_1000g_pon.hg38.vcf.gz
resource=$root/data/known/somatic-hg38_af-only-gnomad.hg38.vcf.gz
res_dir=$root/result/ref_mutect2/$sample
interval=$root/data/known/$2.interval_list

# use origin and restored bams
truncated=$root/script/callSV/truncated_bams
origin="/lustre/home/acct-clswcc/clswcc-zqd/project/GATK_bam"
restore=$root/data/gatkBams
if [[ $(cat $truncated) =~ $normal ]]; then
	normal_bam=$restore/$normal/$normal.recal.bam
else
	normal_bam=$origin/$normal/$normal.recal_reads.bam
fi

if [[ $(cat $truncated) =~ $tumor ]]; then
	tumor_bam=$restore/$tumor/$tumor.recal.bam
else
	tumor_bam=$origin/$tumor/$tumor.recal_reads.bam
fi

# GATK mutect2 pipeline
gatk Mutect2 \
  -R $reference \
  -I $tumor_bam \
  -I $normal_bam \
  -tumor $tumor \
  -normal $normal \
  -L $interval \
  -pon $pon \
  --germline-resource $resource \
  -O $res_dir/${sample}_${2}.vcf.gz

