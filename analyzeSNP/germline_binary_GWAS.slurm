#!/bin/bash

#SBATCH -J germline_GWAS_case_control
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/snp
res_dir=$root/result/snpGWAS/binary

cd $res_dir || exit

control=("hb" "hs")
for c in ${control[*]} ; do
  # Convert vcf file to plink bed/bim/fam files
  plink --vcf $vcf_dir/germline_${c}_snp.vcf --out germline_${c}_binary --double-id

  # Add phenotype info
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f germline_${c}_binary.fam \
    -o germline_${c}_binary.fam \
    -t "Case"
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/hanBeijing1000genomes/${c}_pheno.txt \
    -f germline_${c}_binary.fam \
    -o germline_${c}_binary.fam \
    -t "Control"

  # Plink QC
  # Step 1 Delete SNPs and individuals with high levels of missingness
  plink --bfile germline_${c}_binary --geno 0.02 --make-bed --out germline_tmp_1
  plink --bfile germline_tmp_1 --mind 0.02 --make-bed --out germline_tmp_2
  # Step 2 Remove SNPs with a low MAF frequency
  plink --bfile germline_tmp_2 --maf 0.01 --make-bed --out germline_tmp_3
  # Step 3 Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
  plink --bfile germline_tmp_3 --hwe 1e-6 --make-bed --out germline_tmp_4
  plink --bfile germline_tmp_4 --hwe 1e-6 --make-bed --out germline_tmp_5
  plink --bfile germline_tmp_5 --hwe 1e-10 --hwe-all --make-bed --out germline_${c}_binary_qc

  # Plink assoc
  plink --bfile germline_${c}_binary_qc -assoc --adjust --out germline_${c}_binary_assoc
  # Plink logistic
  plink --bfile germline_${c}_binary_qc --logistic --out germline_${c}_binary_logistic
  awk '!/"NA"/' germline_${c}_binary_logistic.assoc.logistic > germline_${c}_binary_logistic_no_default.assoc.logistic
done
