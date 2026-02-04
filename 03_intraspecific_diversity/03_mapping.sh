#!/bin/bash
# 07-10a_map_picard_samptools.slurm
# This script submits a batch of read-mapping processes using SLURM.

# Loop through each sample from the list
for sample in $(cut -f1 /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/ahya_gr_sample.list | grep 'RRAP-ECO3-2021-Ahya-LIZA-446_L4'); do
    sbatch --partition=cpuq \
    --job-name=mapAhGR_${sample} \
    --output=/export/home/a-e/amateos/phyluce_test_2/mapping_tests/logs/ahya_GR_map_logs/missWG_${sample}.log \
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/aspat_bowtie2_picard_samtools.slurm \
        $sample \
        /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/02-05_clean_reads \
        /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/ahya_gr_meta.list \
        /export/home/a-e/amateos/phyluce_test_2/mapping_tests/reference/bowtie_databases \
        'ahya_WG' \
        /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/07_map_reads/
done

# $1 = sample, $2 = indir, $3 = infolist, $4 = refdir, $5 = reflist, $6 = outdir

##########| grep 'x' for single sample test in line 6

### RRAP-ECT01-2021-Ahya-B6_L3
### RRAP-ECO3-2021-Ahya-MOOR-428_L3
### RRAP-ECO3-2021-Ahya-LIZA-446_L4



