#!/bin/bash

#SBATCH -J indel_germline_binary
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/indel
res_dir=$root/result/indelGWAS/binary
sum_dir=$root/result/indelGWAS/summary

cd $res_dir || exit

control=("hb" "hs")
for c in ${control[*]} ; do
  # Convert vcf file to plink bed/bim/fam files
  plink --vcf $vcf_dir/germline_${c}.vcf --out germline_${c}_binary --double-id

  # Add phenotype info
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f $res_dir/germline_${c}_binary.fam \
    -o $res_dir/germline_${c}_binary.fam \
    -t "Case"
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/${c}_pheno.txt \
    -f $res_dir/germline_${c}_binary.fam \
    -o $res_dir/germline_${c}_binary.fam \
    -t "Control"

  # Plink QC
  # Step 1 Delete SNPs and individuals with high levels of missingness
  plink --bfile germline_${c}_binary --geno 0.02 --make-bed --out germline_${c}_tmp_1
  plink --bfile germline_${c}_tmp_1 --mind 0.02 --make-bed --out germline_${c}_tmp_2
  # Step 2 Remove SNPs with a low MAF frequency
  plink --bfile germline_${c}_tmp_2 --maf 0.01 --make-bed --out germline_${c}_tmp_3
  # Step 3 Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
  plink --bfile germline_${c}_tmp_3 --hwe 1e-6 --make-bed --out germline_${c}_tmp_4
  plink --bfile germline_${c}_tmp_4 --hwe 1e-10 --hwe-all --make-bed --out germline_${c}_binary_qc

  # Plink assoc
  plink --bfile germline_${c}_binary_qc -assoc --adjust --out germline_${c}_binary_assoc
  # Plink logistic
  plink --bfile germline_${c}_binary_qc --logistic --adjust --out germline_${c}_binary_logistic
  awk '!/ NA/' germline_${c}_binary_logistic.assoc.logistic > ${sum_dir}/${c}

  # Delete temporary files
  rm germline_${c}_tmp_*.b*
  rm germline_${c}_tmp_*.f*
  rm germline_${c}_tmp_*.h*
done

phenotype=("EBV" "HP")
for p in ${phenotype[*]} ; do
  mkdir $res_dir/$p
  cd $res_dir/$p || exit

  plink --vcf $vcf_dir/germline.vcf --out germline_$p --double-id
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f $res_dir/$p/germline_$p.fam \
    -o $res_dir/$p/germline_$p.fam \
    -t $p
  sed -i "s/1$/2/g" germline_$p.fam
  sed -i "s/0$/1/g" germline_$p.fam
  plink --bfile germline_$p --geno 0.02 --make-bed --out germline_${p}_tmp_1
  plink --bfile germline_${p}_tmp_1 --mind 0.02 --make-bed --out germline_${p}_tmp_2
  plink --bfile germline_${p}_tmp_2 --maf 0.01 --make-bed --out germline_${p}_tmp_3
  plink --bfile germline_${p}_tmp_3 --hwe 1e-6 --make-bed --out germline_${p}_tmp_4
  plink --bfile germline_${p}_tmp_4 --hwe 1e-10 --hwe-all --make-bed --out germline_${p}_qc

  plink --bfile germline_${p}_qc -assoc --adjust --out germline_${p}_assoc
  plink --bfile germline_${p}_qc --logistic --out germline_${p}_logistic
  awk '!/NA/' germline_${p}_logistic.assoc.logistic > ${sum_dir}/${p}
done