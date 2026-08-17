#!/bin/bash

# 06_fastqc_post.sh
# This script submits one FastQC job per sample.
# After this script is complete, run the following commands to get a MultiQC summary:
#   cd result/06_fastqc
#   multiqc *


for sample in $(cut -f1 /your/path/ahya_sample.list); do
    sbatch --job-name=fastqc_${sample} \
           --output=/your/path/fqc_${sample}.out \
           /your/path/fastqc.slurm \
           "/your/path/results/02-05_clean_reads/${sample}_SingleEndQualFiltered_norm.fastq.gz" \
           "/your/path/results/02-05_clean_reads/${sample}.notCombined_norm_1.fastq.gz" \
           "/your/path/results/02-05_clean_reads/${sample}.notCombined_norm_2.fastq.gz" \
           "/your/path/results/06_fastqc/"
done
