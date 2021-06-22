#!/bin/bash

#SBATCH -J indel_germline_categorical
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH --output=%J.out
#SBATCH --error=%J.err

root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
spt_dir=$root/script/analyzeIndel
vcf_dir=$root/data/indel
res_dir=$root/result/indelGWAS/categorical

cd $res_dir || exit

# Convert vcf file to plink bed/bim/fam files
plink --vcf $vcf_dir/germline.vcf --out germline --double-id

# Add gender information
python3 $spt_dir/add_phenotype.py \
  -p $vcf_dir/phenotype.txt \
  -f $res_dir/germline.fam \
  -o $res_dir/germline.fam \
  -t "Location"

# Plink QC
# Step 1 Delete SNPs and individuals with high levels of missingness
plink --bfile germline --geno 0.02 --make-bed --out germline_tmp_1
plink --bfile germline_tmp_1 --mind 0.02 --make-bed --out germline_tmp_2
# Step 2 Remove SNPs with a low MAF frequency
plink --bfile germline_tmp_2 --maf 0.01 --make-bed --out germline_tmp_3
# Step 3 Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE)
plink --bfile germline_tmp_3 --hwe 1e-10 --hwe-all --make-bed --out germline_qc

# Delete temporary files
rm germline_tmp_*.b*
rm germline_tmp_*.f*
rm germline_tmp_*.h*

# Run GWAS for different phenotypes
phenotypes=("Subtype" "Location" "Borrman" "Grade" "Lauren")
for p in "${phenotypes[@]}" ; do
  # make phenotype directory
  mkdir $res_dir/$p
  cp $res_dir/germline_qc.b* $res_dir/$p/

  # Add phenotype info to fam file
  python3 $spt_dir/add_phenotype.py \
    -p $vcf_dir/phenotype.txt \
    -f $res_dir/germline_qc.fam \
    -o $res_dir/$p/germline_qc.fam \
    -s " " \
    -t $p

  # GWAS for categorical phenotypes

done

