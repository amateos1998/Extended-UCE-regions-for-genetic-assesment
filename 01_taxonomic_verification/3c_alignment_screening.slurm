#!/bin/bash
#SBATCH --job-name=mafftall    # Job name
#SBATCH --partition=cpuq   # Specify the partition you want to use
#SBATCH --cpus-per-task=20                 # Number of CPU cores per task
#SBATCH --time=48:00:00                   # Time limit hrs:min:sec
#SBATCH --mem=100G                          # Memory per node


# Clean  mafft alignments (remove locus names)

phyluce_align_remove_locus_name_from_files \
    --alignments $INPUT_DIR \
    --output $OUTPUT_DIR \
    --input-format fasta \
    --output-format nexus \
    --cores 8 \
    --log-path logs

 

# # --- Alignment Screening ---
# # Screen alignments for problems before any trimming
# echo "--- Screening Alignments for Problems ---"

phyluce_align_screen_alignments_for_problems \
    --alignments aln-clean \
    --output screened-alignments \
    --input-format nexus \
    --log-path logs \
    --cores 20


# # --- Remove Empty Taxa ---

# # Remove taxa with no data from the screened alignments

# echo "--- Removing Empty Taxa ---"

phyluce_align_remove_empty_taxa \
    --alignments screened-alignments \
    --output screened-no-empty-taxa \
    --input-format nexus \
    --output-format nexus \
    --log-path logs \
    --cores 20