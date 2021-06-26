#!/bin/bash

#SBATCH -J pav_somatic_quantitative
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/PAV
res_dir=$root/result/pavGWAS/quantitative

cd $res_dir || exit

# Copy plink files after qc from categorical phenotypes GWAS
cp ../categorical/somatic_qc.* ./

# Add phenotype information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/somatic_qc.fam \
  -o $res_dir/somatic_qc.fam \
  -t "Diameter"

# Plink assoc
plink --bfile somatic_qc -assoc --adjust --allow-extra-chr --out somatic_assoc
# Plink linear
plink --bfile somatic_qc --linear --allow-extra-chr --out somatic_linear
awk '!/NA/' somatic_linear.assoc.linear > somatic_linear_no_default.assoc.linear