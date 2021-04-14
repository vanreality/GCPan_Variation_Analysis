#!/bin/bash

#SBATCH -J blastn
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=8

fasta=$1
result=$2
db="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/novelSeq/pan_blastn_idx"
blastn="/lustre/home/acct-clswcc/clswcc-fsy/software/ncbi-blast-2.10.1+/bin/blastn"

$blastn -db $db \
  -query $fasta \
  -out $result \
  -outfmt 6 \
  -num_threads 8 \
  -evalue 1e-5 \
  -max_target_seqs 100