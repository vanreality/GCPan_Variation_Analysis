#!/bin/bash

#SBATCH -J somatic_GWAS
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeSNP
vcf_dir=$root/data/snp
res_dir=$root/result/snpGWAS/categorical

cd $res_dir || exit

# Convert vcf file to plink bed/bim/fam files
plink --vcf $vcf_dir/somatic_mutect2_snp.vcf --allow-extra-chr --out somatic_mutect2 --double-id

# Add gender information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/somatic_mutect2.fam \
  -o $res_dir/somatic_mutect2.fam \
  -t "Location"

# Plink QC
# Step 1 Delete SNPs and individuals with high levels of missingness
plink --bfile somatic_mutect2 --geno 0.02  --allow-extra-chr --make-bed --out somatic_tmp_1
plink --bfile somatic_tmp_1 --mind 0.02  --allow-extra-chr --make-bed --out somatic_tmp_2
# Step 2 Remove SNPs with a low MAF frequency
plink --bfile somatic_tmp_2 --maf 0.01  --allow-extra-chr --make-bed --out somatic_tmp_3
# Step 3 Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
plink --bfile somatic_tmp_3 --hwe 1e-10  --allow-extra-chr --hwe-all --make-bed --out somatic_mutect2_qc

# Delete temporary files
rm somatic_tmp_*.b*
rm somatic_tmp_*.f*

# Run GWAS for different phenotypes
phenotypes=("Location" "Borrman" "Grade" "Lauren" "EBV" "HP")
for p in "${phenotypes[@]}" ; do
  # make phenotype directory
  mkdir $res_dir/$p
  cp $res_dir/somatic_mutect2_qc.b* $res_dir/$p/

  # Add phenotype info to fam file
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f $res_dir/somatic_mutect2_qc.fam \
    -o $res_dir/$p/somatic_mutect2_qc.fam \
    -s " " \
    -t $p

  # GWAS for categorical phenotypes

done