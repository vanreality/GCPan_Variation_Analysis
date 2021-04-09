#!/bin/bash

#SBATCH -J svaba_convert
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=convert.out
#SBATCH --error=convert.err
#SBATCH --ntasks-per-node=1

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
svaba=$root/result/svaba
convert=$root/script/callSV/svaba_convert.R

samples=$(cat 185samples)
for sample in $samples ; do
    cd $svaba/$sample || exit
    singularity run /lustre/share/img/ood-app/rstudio_3.6.3.sif Rscript $convert $svaba/$sample/$sample.svaba.germline.sv.vcf
    singularity run /lustre/share/img/ood-app/rstudio_3.6.3.sif Rscript $convert $svaba/$sample/$sample.svaba.somatic.sv.vcf
done
