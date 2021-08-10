#!/bin/bash

#SBATCH -J ref_fiter
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
known_dir=$root/data/known
chrom_list=$(cat $root/script/callSNP/chrom_list)

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

# GATK merge somatic variants of each chromosome
zcat ${sample}_chr1.vcf.gz | grep "^#" > ${sample}_merged.vcf
for i in $chrom_list ; do
  zcat ${sample}_${i}.vcf.gz | grep -v "^#" >> ${sample}_merged.vcf
done

gzip ${sample}_merged.vcf && rm ${sample}_merged.vcf

# GATK GetPileupSummaries to summarize read support for a set number of known variant sites
gatk GetPileupSummaries \
  -I $tumor_bam \
  -V $known_dir/chr17_small_exac_common_3_grch38.vcf.gz \
  -L $known_dir/chr17_small_exac_common_3_grch38.vcf.gz \
  -O tumor_get_pileup_summaries.table

gatk GetPileupSummaries \
  -I $normal_bam \
  -V $known_dir/chr17_small_exac_common_3_grch38.vcf.gz \
  -L $known_dir/chr17_small_exac_common_3_grch38.vcf.gz \
  -O normal_get_pileup_summaries.table

# Estimate contamination
gatk CalculateContamination \
  -I tumor_get_pileup_summaries.table \
  -matched normal_get_pileup_summaries.table \
  -O calculate_contamination.table

# Apply filters
gatk FilterMutectCalls \
  -V ${sample}_merged.vcf.gz \
  --contamination-table calculate_contamination.table \
  --stats filtered.stats.txt \
  -O ${sample}_filtered.vcf.gz