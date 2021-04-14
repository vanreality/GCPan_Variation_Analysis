#!/bin/bash

#SBATCH -J insertion
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=1

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzeSV
manta_dir=$root/result/manta

samples=$(cat $script_dir/185samples)
for sample in $samples; do
  vcf_dir=$manta_dir/$sample/result/results/variants
  singularity run /lustre/share/img/ood-app/rstudio_3.6.3.sif Rscript $script_dir/manta_extract_insertions.R $vcf_dir/germlineSV.vcf
  singularity run /lustre/share/img/ood-app/rstudio_3.6.3.sif Rscript $script_dir/manta_extract_insertions.R $vcf_dir/somaticSV.vcf
done;