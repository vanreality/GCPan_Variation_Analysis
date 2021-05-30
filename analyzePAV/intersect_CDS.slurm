#!/bin/bash

#SBATCH -J intersect
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzePAV
absence_region=$root/result/PAV/absenceRegion
intersect=$root/result/PAV/intersect
cds=$root/data/PAV/cds.ann.bed
vcf_dir=$root/result/survivor

cat < $script_dir/185samples | while IFS= read -r sample
do
  normal=$(echo $sample | awk -F "." '{print $1}')
  tumor=$(echo $sample | awk -F "." '{print $2}')

  pav_normal=$absence_region/$normal.absence.bed
  pav_tumor=$absence_region/$tumor.absence.bed

  sv_normal=$vcf_dir/$sample/$sample.germline.merged.vcf.bed
  sv_tumor=$vcf_dir/$sample/$sample.somatic.merged.vcf.bed

  bedtools intersect -a $sv_normal -b $pav_normal -wa -wb > $intersect/$normal.bed
  bedtools intersect -a $intersect/$normal.bed -b $cds -wa -wb > $intersect/$normal.cds.bed

  bedtools intersect -a $sv_normal -b $pav_tumor -wa -wb > $intersect/$tumor.bed
  bedtools intersect -a $sv_tumor -b $pav_tumor -wa -wb >> $intersect/$tumor.bed
  bedtools intersect -a $intersect/$tumor.bed -b $cds -wa -wb > $intersect/$tumor.cds.bed
done
