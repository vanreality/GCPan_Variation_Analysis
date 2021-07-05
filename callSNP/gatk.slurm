#!/bin/bash

#SBATCH -J gatk
#SBATCH -p arm128c256g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

# arm software
module use /lustre/share/spack/modules/kunpeng920/linux-centos7-aarch64
module load gatk/4.2.0.0-gcc-9.3.0-openblas
module load bwa/0.7.17-gcc-9.3.0
samtools=/lustre/home/acct-clswcc/clswcc-xhz/armtools/samtools-1.12/src/bin/samtools

# fastq, reference, known snp/indel
sample=$1
reference="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/known/pan.fa"
fastq_file1="/lustre/home/acct-clswcc/clswcc-zqd/project/fastq/"${sample}_combined_R1.fastq.gz
fastq_file2="/lustre/home/acct-clswcc/clswcc-zqd/project/fastq/"${sample}_combined_R2.fastq.gz
knownsites_dir="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/known"

#step 1. Mapping 'fastq' file to the reference by bwa:
echo -e "Job $SLURM_JOBID starts at $(date)." >>gatk.log
start=$(date +%s)
bwa mem -M -R "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:bar\tPU:barcode" \
	-t 16 $reference $fastq_file1 $fastq_file2 >$sample.sam 2>$sample.log
end=$(date +%s)
cost=$(( end - start ))
echo -e "bwa kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log

start=$(date +%s)
$samtools view -@ 16 -bS $sample.sam > $sample.bam \
   && rm $sample.sam
end=$(date +%s)
cost=$(( end - start ))
echo -e "sam2bam kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log

start=$(date +%s)
$samtools sort -@ 16 -m 32G -o $sample.sorted.bam -O bam -T sorted $sample.bam \
   && rm $sample.bam
end=$(date +%s)
cost=$(( end - start ))
echo -e "sortbam kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log

#step 2. Mark PCR duplications by gatk:
start=$(date +%s)
gatk MarkDuplicates \
	-I $sample.sorted.bam \
	-O $sample.dedup.bam \
	-M metrics.txt \
   && rm $sample.sorted.bam
end=$(date +%s)
cost=$(( end - start ))
echo -e "MarkDuplicates kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log

start=$(date +%s)
$samtools index -@ 16 $sample.dedup.bam
end=$(date +%s)
cost=$(( end - start ))
echo -e "index kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log

#step 3. Base quality score recalibration by gatk:
start=$(date +%s)
gatk BaseRecalibrator \
   -R $reference \
   -I $sample.dedup.bam \
   --known-sites $knownsites_dir/dbsnp_146.hg38.vcf \
   --known-sites $knownsites_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   -O recal_data.table

gatk ApplyBQSR \
   -R $reference \
   -I $sample.dedup.bam \
   -bqsr recal_data.table \
   -O $sample.recal.bam \
   && rm $sample.dedup.ba*
end=$(date +%s)
cost=$(( end - start ))
echo -e "BQSR kernel time $((cost/60))min $((cost%60))s.\n" >>gatk.log
