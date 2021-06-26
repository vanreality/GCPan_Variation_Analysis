#!/bin/bash

#SBATCH -J pav_somatic_categorical
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/PAV
res_dir=$root/result/pavGWAS/categorical

cd $res_dir || exit

# Convert vcf file to plink bed/bim/fam files
plink --vcf $vcf_dir/tumor.vcf --out somatic --double-id --allow-extra-chr

# Add gender information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/somatic.fam \
  -o $res_dir/somatic.fam \
  -t "Location"

# Plink QC
# Step 1 Delete SNPs and individuals with high levels of missingness
plink --bfile somatic --geno 0.02 --allow-extra-chr --make-bed --out somatic_tmp_1
plink --bfile somatic_tmp_1 --mind 0.02 --allow-extra-chr --make-bed --out somatic_tmp_2
# Step 2 Remove SNPs with a low MAF frequency
plink --bfile somatic_tmp_2 --maf 0.01 --allow-extra-chr --make-bed --out somatic_tmp_3
# Step 3 Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
plink --bfile somatic_tmp_3 --hwe 1e-10 --allow-extra-chr --hwe-all --make-bed --out somatic_qc

# Delete temporary files
rm somatic_tmp_*.b*
rm somatic_tmp_*.f*
rm somatic_tmp_*.h*

# Run GWAS for different phenotypes
phenotypes=("Subtype" "Location" "Grade" "Lauren")
for p in ${phenotypes[*]} ; do
  # make phenotype directory
  mkdir $res_dir/$p
  cp $res_dir/somatic_qc.b* $res_dir/$p/
  cd $res_dir/$p || exit

  # Add phenotype info to fam file
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f $res_dir/somatic_qc.fam \
    -o $res_dir/$p/somatic_qc.fam \
    -s " " \
    -t $p

  # GWAS for categorical phenotypes
  categories=$(awk '{print $6}' somatic_qc.fam | sort | uniq)
  for c in ${categories[*]}; do
    mkdir $res_dir/$p/$c
    cd $res_dir/$p/$c || exit
    cp $res_dir/$p/somatic_qc* ./
    bash $spt_dir/case_control.sh somatic_qc.fam $c

    # Treat every subtype as binary trait
    # Plink assoc
    plink --bfile somatic_qc -assoc --adjust --allow-extra-chr --out somatic_assoc
    # Plink logistic
    plink --bfile somatic_qc --logistic --allow-extra-chr --out somatic_logistic
    awk '!/NA/' somatic_logistic.assoc.logistic > somatic_logistic_no_na.assoc.logistic
  done
done

