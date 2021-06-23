#!/bin/bash

#SBATCH -J indel_germline_quantitative
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeIndel
vcf_dir=$root/data/indel
res_dir=$root/result/indelGWAS/quantitative

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
plink --bfile germline_qc --linear --out germline_linear
awk '!/NA/' germline_linear.assoc.linear > germline_linear_no_default.assoc.linear