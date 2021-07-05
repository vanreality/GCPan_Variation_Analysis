#!/bin/bash

#SBATCH -J snp_germline_quantitative
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/snp
res_dir=$root/result/snpGWAS/quantitative
sum_dir=$root/result/snpGWAS/summary

cd $res_dir || exit

# Copy plink files after qc from categorical phenotypes GWAS
cp ../categorical/germline_qc.* ./

# Add phenotype information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/germline_qc.fam \
  -o $res_dir/germline_qc.fam \
  -t "Diameter"

# Plink assoc
plink --bfile germline_qc -assoc --adjust --out germline_assoc
# Plink linear
plink --bfile germline_qc --linear --adjust --out germline_linear
awk '!/NA/' germline_linear.assoc.linear > ${sum_dir}/Diameter