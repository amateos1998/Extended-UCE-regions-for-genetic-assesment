#!/bin/bash
#SBATCH --job-name=ed_tr
#SBATCH --output=edgetrim_%j.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=9-24:00:00
#SBATCH --partition=cpuq


date
echo
echo "edge trimming alignments"
echo

phyluce_align_get_trimmed_alignments_from_untrimmed \
    --alignments screened2 \
    --output edge-trimmed \
    --input-format fasta \
    --cores 20

echo
echo "edge trimming complete"
date

echo "get summary data for edge trimmed aligs"
phyluce_align_get_align_summary_data \
    --alignments edge-trimmed \
    --cores 20

 
