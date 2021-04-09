#!/bin/bash

#SBATCH -J survivor
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=survivor.out
#SBATCH --error=survivor.err
#SBATCH --ntasks-per-node=1

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
survivor="/lustre/home/acct-clswcc/clswcc-fsy/software/SURVIVOR-master/Debug/SURVIVOR"
result=$root/result/survivor
delly=$root/result/delly
svaba=$root/result/svaba
manta=$root/result/manta
samples=$(cat 185samples)
#samples="WGC051030D.WGC051031D"

if [ ! -d $result ]; then
    mkdir $result
fi
cd $result || exit

for sample in $samples ; do
  if [ -d $sample ]; then
      continue
  fi
  mkdir $sample
  cd $sample || exit

  # merge germline vcf files
  {
    printf "%s/%s/result/results/variants/germlineSV.vcf\n" $manta $sample
    printf "%s/%s/%s.vcf\n" $delly $sample $sample
    printf "%s/%s/%s.svaba.germline.sv.vcf.converted" $svaba $sample $sample
  } > germline_vcf_files

  # merge somatic vcf files
  {
    printf "%s/%s/result/results/variants/somaticSV.vcf\n" $manta $sample
    printf "%s/%s/%s.somatic.vcf\n" $delly $sample $sample
    printf "%s/%s/%s.svaba.somatic.sv.vcf.converted" $svaba $sample $sample
  } > somatic_vcf_files

  $survivor merge germline_vcf_files 1000 2 1 1 0 30 ${sample}.germline.merged.vcf

  $survivor merge somatic_vcf_files 1000 2 1 1 0 30 ${sample}.somatic.merged.vcf

  cd ..
done