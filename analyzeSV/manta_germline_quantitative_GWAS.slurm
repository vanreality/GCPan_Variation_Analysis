#!/bin/bash

#SBATCH -J sv_germline_manta_quantitative
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/result/survivor/survivor_germline
res_dir=$root/result/svGWAS/quantitative

cd $res_dir || exit

# Copy plink files after qc from categorical phenotypes GWAS
cp ../categorical/germline_manta_qc.* ./

# Add phenotype information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/germline_manta_qc.fam \
  -o $res_dir/germline_manta_qc.fam \
  -t "Diameter"

# Plink assoc
plink --bfile germline_manta_qc -assoc --allow-extra-chr --adjust --out germline_manta_assoc
# Plink linear
plink --bfile germline_manta_qc --linear --allow-extra-chr --out germline_manta_linear
awk '!/NA/' germline_manta_linear.assoc.linear > germline_manta_linear_no_default.assoc.linear