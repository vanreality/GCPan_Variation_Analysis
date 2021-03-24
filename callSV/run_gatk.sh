#!/bin/bash

script="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/callSV"
result="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/gatkBams"

samples=`cat truncated_bams`
mkdir $result
cd $result
for sample in $samples;
do
if [ -d $sample ]; then
continue
fi
mkdir $sample;
cd $sample;
cp $script/gatk.slurm ./

sbatch gatk.slurm $sample

cd ..;
done;
