#!/bin/bash

#01_seqprep.sh
#This script subimts jobs to prepare raw reads using bin/bbsplit_clean_norm.slurm
#All samples will be performed for bbsplit size-normalised read-QC and decontamination
#You will need a .list file with your sample names + metadata

for sample in $(cut -f1 /your/path/sample.list); do
  sbatch --partition=cpuq \
  --job-name=seqprep_${sample} \
  --output=/your/path/_${sample}.out \
  /your/path/bbsplit_clean_norm.slurm $sample your/path/01_results
done




