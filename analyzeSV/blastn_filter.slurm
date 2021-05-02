#!/bin/bash

#SBATCH -J blastn_filter
#SBATCH -p small
#SBATCH --mail-user=mr.van@sjtu.edu.cn
#SBATCH --mail-type=ALL
#SBATCH --output=%J.out
#SBATCH --error=%J.err
#SBATCH --ntasks-per-node=2

python3 blastn_filter.py -i $1 -f $2