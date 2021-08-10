#!/bin/bash

#SBATCH -J manta
#SBATCH -p arm128c256g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=16

module use /lustre/share/spack/modules/kunpeng920/linux-centos7-aarch64
module load manta/1.6.0-gcc-9.3.0

normal=$1
tumor=$2
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
reference=$root/data/known/pan.fa
normal_bam=$root/data/panGatkBams/$normal/$normal.recal.bam
tumor_bam=$root/data/panGatkBams/$tumor/$tumor.recal.bam

configManta.py --normalBam $normal_bam --tumorBam $tumor_bam --referenceFasta $reference --runDir result
$PWD/result/runWorkflow.py -m local -j 16
