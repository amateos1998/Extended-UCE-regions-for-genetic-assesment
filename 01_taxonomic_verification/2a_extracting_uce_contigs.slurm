#!/bin/bash
#SBATCH --job-name=aHYA_ex     # Job name
#SBATCH --partition=cpuq   # Specify the partition you want to use
#SBATCH --cpus-per-task=20                 # Number of CPU cores per task
#SBATCH --time=48:00:00                   # Time limit hrs:min:sec
#SBATCH --mem=200G                          # Memory per node
#SBATCH --output=/export/home/a-e/amateos/phyluce_test_2/logs/extractUCE_%j.out  # Standard output log

###THIS SCRIPT extracts UCE loci from contigs assembled from lcWGS

# Activate the conda environment where Phyluce is installed
source /export/home/a-e/amateos/miniconda3/etc/profile.d/conda.sh
conda activate phyluce-1.7.3



date
echo "get match counts"

phyluce_assembly_get_match_counts \
    --locus-db /export/home/a-e/amateos/phyluce_test_2/results/apples_ahya_uces/uce_search_results/probe.matches.sqlite \
    --taxon-list-config /export/home/a-e/amateos/phyluce_test_2/results/apples_ahya_uces/uce_search_results/taxa.conf \
    --taxon-group 'all' \
    --incomplete-matrix \
    --output all-taxa-incomplete.conf
echo
date
done