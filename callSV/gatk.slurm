#!/bin/bash

#SBATCH -J gatk
#SBATCH -p small
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=8

#module load bwa
module load samtools

sample=$1
reference="/lustre/home/acct-clswcc/clswcc-zqd/data/hg38/hg38.fa"
fastq_file="/lustre/home/acct-clswcc/clswcc-zqd/project/fastq/"$sample*
knownsites_dir="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/known"
bwa="/lustre/home/acct-clswcc/clswcc-zqd/software/bwa/bwa-0.7.10/bwa"
gatk="/lustre/home/acct-clswcc/clswcc-fsy/software/gatk-4.1.9.0/gatk"

#step 1. Mapping 'fastq' file to the reference by bwa:
echo -e "Job $SLURM_JOBID starts at $(date)." >>gatk.log
$bwa mem -M -R "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:bar\tPU:barcode" \
	-t 8 $reference $fastq_file >$sample.sam 2>$sample.log
echo -e "bwa end at $(date).\n" >>gatk.log

echo -e "sam2bam start at $(date)." >>gatk.log
samtools view -@ 8 -bS $sample.sam > $sample.bam \
   && rm $sample.sam
echo -e "sam2bam end at $(date).\n" >>gatk.log

echo -e "sortbam start at $(date)." >>gatk.log
samtools sort -@ 8 -m 4G -o $sample.sorted.bam -O bam -T sorted $sample.bam \
   && rm $sample.bam
echo -e "sortbam end at $(date).\n" >>gatk.log

#step 2. Mark PCR duplications by gatk:
echo -e "MarkDuplicates start at $(date)." >>gatk.log
$gatk MarkDuplicates \
	-I $sample.sorted.bam \
	-O $sample.dedup.bam \
	-M metrics.txt \
   && rm $sample.sorted.bam
echo -e "MarkDuplicates end at $(date).\n" >>gatk.log

echo -e "Index start at $(date)." >>gatk.log
samtools index -@ 8 $sample.dedup.bam
echo -e "Index end at $(date).\n" >>gatk.log

#step 3. Base quality score recalibration by gatk:
echo -e "BQSR start at $(date)." >>gatk.log
$gatk BaseRecalibrator \
   -I $sample.dedup.bam \
   -R $reference \
   --known-sites $knownsites_dir/dbsnp_146.hg38.vcf \
   --known-sites $knownsites_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   -O recal_data.table

$gatk ApplyBQSR \
   -R $reference \
   -I $sample.dedup.bam \
   -bqsr recal_data.table \
   -O $sample.recal.bam \
   && rm $sample.dedup.ba*
echo -e "BQSR end at $(date).\n" >>gatk.log
