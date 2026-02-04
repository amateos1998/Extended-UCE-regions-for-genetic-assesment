#!/bin/bash
#SBATCH --job-name=catphy
#SBATCH --output=catphy_%j.log
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=9-24:00:00
#SBATCH --partition=cpuq

echo
date

echo "# Get concatenated matrices for edge-trimmed alignments
2

#    a. 50% complete matrix

phyluce_align_get_only_loci_with_min_taxa \
    --alignments edge-trimmed \
    --taxa 963 \ #no. of loci inside edge-trimmed folder
    --percent 0.50 \  #% of completeness of matrix (50, 75, 90 etc.)
    --output edge-trimmed-50p \
    --cores 8 \
    --log-path logs

 

#    b. 75% complete matrix

phyluce_align_get_only_loci_with_min_taxa \
    --alignments edge-trimmed \
    --taxa 963 \
    --percent 0.75 \
    --output edge-trimmed-75p \
    --cores 8 \
    --log-path logs

echo
done
date


#optional to get very large alig
#phyluce_align_concatenate_alignments \
 #   --alignments \
#--input-format fasta \
 #   --output \
  #  --phylip \
   # --log-path log
