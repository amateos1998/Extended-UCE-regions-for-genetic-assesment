#!/bin/bash
# 07-10a_map_picard_samptools.slurm
# This script submits a batch of read-mapping processes using SLURM.

# Loop through each sample from the list
for sample in $(cut -f1/your/path/ahya_gr_sample.list | grep 'RRAP-ECO3-2021-Ahya-LIZA-446_L4'); do
    sbatch --partition=cpuq \
    --job-name=mapAhGR_${sample} \
    --output=/your/path/_${sample}.log \
   /your/path/bowtie2_picard_samtools.slurm \
        $sample \
        /your/path/results/02-05_clean_reads \
       /your/path/data/ahya_gr_meta.list \
        /your/path/reference/bowtie_databases \
        'ahya_WG' \
       /your/path/results/07_map_reads/
done




