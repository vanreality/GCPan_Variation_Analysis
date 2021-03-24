#!/bin/bash

script="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSV"
result="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/svaba"

samples=`cat 185samples`
mkdir $result
cd $result
for sample in $samples;
do
if [ -d $sample ]; then
continue
fi
mkdir $sample;
cd $sample;

#batch process

normal=$(echo $sample | awk -F "." '{print $1}')
tumor=$(echo $sample | awk -F "." '{print $2}')

cp $script/svaba.slurm ./
sbatch svaba.slurm $normal $tumor

cd ..;
done;
