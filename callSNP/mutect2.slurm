#!/bin/bash

#SBATCH -J gatk
#SBATCH -p arm128c256g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=1

# arm software
module use /lustre/share/spack/modules/kunpeng920/linux-centos7-aarch64
module load gatk/4.2.0.0-gcc-9.3.0-openblas
module load bwa/0.7.17-gcc-9.3.0
module load gcc
samtools=/lustre/home/acct-clswcc/clswcc-xhz/armtools/samtools-1.12/src/bin/samtools

# data
sample=$1
normal=$(echo $sample | awk -F "." '{print $1}')
tumor=$(echo $sample | awk -F "." '{print $2}')
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
reference=$root/data/known/pan.fa
normal_bam=$root/data/panGatkBams/$normal/$normal.recal.bam
tumor_bam=$root/data/panGatkBams/$tumor/$tumor.recal.bam
pon=$root/data/known/somatic-hg38_1000g_pon.hg38.vcf.gz
resource=$root/data/known/somatic-hg38_af-only-gnomad.hg38.vcf.gz
res_dir=$root/result/mutect2/$sample
interval=$root/data/known/$2.interval_list

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

