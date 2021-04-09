#!/bin/bash

#SBATCH -J compare
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=1

module load bedtools2

# path
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script_dir=$root/script/analyzeSV
cds=$root/data/PAV/cds.ann.bed
result=$root/result/compareSVPAV

# python scripts
calculate=$script_dir/calculate_coverage.py
compare=$script_dir/compare_with_PAV.py

#sv_calling_method=("manta" "svaba")
variation_type=("germline" "somatic")

for type in ${variation_type[*]} ; do
  vcf_dir=$root/result/survivor
  res_dir=$result/survivor.$type.sv
  mkdir $res_dir
  mkdir $res_dir/intersect

  cat < $script_dir/185samples | while IFS= read -r sample
  do
    vcf=$vcf_dir/$sample/$sample.${type}.merged.vcf

    singularity run /lustre/share/img/ood-app/rstudio_3.6.3.sif Rscript $script_dir/survivor_extract_gene_sv.R $vcf
    bedtools intersect -a ${vcf}.bed -b $cds -wa -wb > ${res_dir}/intersect/${sample}

    if [[ $type == "somatic" ]]; then
      cat $result/survivor.germline.sv/intersect/${sample} >> $result/survivor.somatic.sv/intersect/${sample}
    fi
  done

  if [[ $type == "germline" ]]; then
    sample_name=$root/data/PAV/normal.individual.info
    PAV=$root/data/PAV/final.normal.cds.cov
  else
    sample_name=$root/data/PAV/tumor.individual.info
    PAV=$root/data/PAV/final.tumor.cds.cov
  fi

  python $calculate -i ${res_dir}/intersect -n $sample_name
  python $compare -s ${res_dir}/coverage -p $PAV
done

# unfiltered

#for method in ${sv_calling_method[*]}; do
#	for type in ${variation_type[*]}; do
#		vcf_dir=$root/result/$method
#		res_dir=$result/$method.$type.sv
#		mkdir $res_dir
#		mkdir $res_dir/intersect
#
#		cat < $script_dir/185samples | while IFS= read -r sample
#		do
#			if [[ $method == "manta" ]]; then
#				vcf=$vcf_dir/$sample/result/results/variants/${type}SV.vcf
#			else
#				vcf=$vcf_dir/$sample/$sample.svaba.${type}.sv.vcf
#			fi
#
#			Rscript $script_dir/${method}_extract_gene_sv.R $vcf
#			bedtools intersect -a ${vcf}.unfiltered.bed -b $cds -wa -wb > ${res_dir}/intersect/${sample}.unfiltered
#
#			if [[ $type == "somatic" ]]; then
#				cat $result/${method}.germline.sv/intersect/${sample}.unfiltered >> $result/${method}.somatic.sv/intersect/${sample}.unfiltered
#			fi
#		done
#
#		if [[ $type == "germline" ]]; then
#			sample_name=$root/data/PAV/normal.individual.info
#			PAV=$root/data/PAV/final.normal.cds.cov
#		else
#			sample_name=$root/data/PAV/tumor.individual.info
#			PAV=$root/data/PAV/final.tumor.cds.cov
#		fi
#
#		python $calculate -i ${res_dir}/intersect -n $sample_name
#		python $compare -s ${res_dir}/coverage.unfiltered -p $PAV
#	done
#done


# filtered with SV length < 100,000

#for method in ${sv_calling_method[*]}; do
#	for type in ${variation_type[*]}; do
#		vcf_dir=$root/result/$method
#		res_dir=$result/$method.$type.sv
#		mkdir $res_dir
#		mkdir $res_dir/intersect
#
#		for sample in $(cat $script_dir/185samples); do
#			if [[ $method == manta ]]; then
#				vcf=$vcf_dir/$sample/result/results/variants/${type}SV.vcf
#			else
#				vcf=$vcf_dir/$sample/$sample.svaba.${type}.sv.vcf
#			fi
#
#			Rscript $script_dir/${method}_extract_gene_sv.R $vcf
#			bedtools intersect -a ${vcf}.bed -b $cds -wa -wb > ${res_dir}/intersect/${sample}
#
#			if [[ $type == somatic ]]; then
#				cat $result/${method}.germline.sv/intersect/${sample} >> $result/${method}.somatic.sv/intersect/${sample}
#			fi
#		done
#
#		if [[ $type == germline ]]; then
#			sample_name=$root/data/PAV/normal.individual.info
#			PAV=$root/data/PAV/final.normal.cds.cov
#		else
#			sample_name=$root/data/PAV/tumor.individual.info
#			PAV=$root/data/PAV/final.tumor.cds.cov
#		fi
#
#		python $calculate -i ${res_dir}/intersect -n $sample_name
#		python $compare -s ${res_dir}/coverage -p $PAV
#	done
#done
