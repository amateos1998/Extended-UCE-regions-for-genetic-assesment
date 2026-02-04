#!/bin/bash

# 06_fastqc_post.sh
# This script submits one FastQC job per sample.
# After this script is complete, run the following commands to get a MultiQC summary:
#   cd result/06_fastqc
#   multiqc *


for sample in $(cut -f1 /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/ahya_sample.list); do
    sbatch --job-name=fastqc_${sample} \
           --output=/export/home/a-e/amateos/phyluce_test_2/mapping_tests/logs/fastqc_all_mod/fqc_${sample}.out \
           /export/home/a-e/amateos/phyluce_test_2/mapping_tests/fastqc.slurm \
           "/export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/02-05_clean_reads/${sample}_SingleEndQualFiltered_norm.fastq.gz" \
           "/export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/02-05_clean_reads/${sample}.notCombined_norm_1.fastq.gz" \
           "/export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/02-05_clean_reads/${sample}.notCombined_norm_2.fastq.gz" \
           "/export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/06_fastqc/"
done