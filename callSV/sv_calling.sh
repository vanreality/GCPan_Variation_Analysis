#!/bin/bash

sv_caller=( manta svaba delly )
if [[ ! (${sv_caller[*]} =~ $1) ]]; then
    printf "Usage: ./run_sv_caller.sh {sv_caller}\nSV caller: manta/svaba/delly\n"
    exit
fi

caller=$1
root="/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
script=$root/script/callSV
result=$root/result/$caller

samples=$(cat 185samples)
if [ ! -d $result ]; then
  mkdir $result
fi
cd $result || exit
for sample in $samples; do
  if [ -d $sample ]; then
    continue
  fi
  mkdir $sample;
  cd $sample || exit

  #batch process

  normal=$(echo $sample | awk -F "." '{print $1}')
  tumor=$(echo $sample | awk -F "." '{print $2}')

  cp $script/$caller.slurm ./
  sbatch $caller.slurm $normal $tumor

  cd ..;
done;
