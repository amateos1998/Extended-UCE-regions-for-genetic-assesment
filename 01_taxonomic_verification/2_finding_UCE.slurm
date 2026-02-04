#!/bin/bash
#SBATCH --job-name=ahyaUCE     # Job name
#SBATCH --partition=cpuq   # Specify the partition you want to use
#SBATCH --cpus-per-task=20                 # Number of CPU cores per task
#SBATCH --time=48:00:00                   # Time limit hrs:min:sec
#SBATCH --mem=200G                          # Memory per node
#SBATCH --output=/export/home/a-e/amateos/phyluce_test_2/logs/findUCE_%j.out  # Standard output log
#SBATCH --error=/export/home/a-e/amateos/phyluce_test_2/logs/findUCE_%j.err 

###THIS SCRIPT FINDS CONTIGS THAT ARE UCE LOCI AND MOVES ASIDE THOSE THAT ARE NOT

# Activate the conda environment where Phyluce is installed
source /export/home/a-e/amateos/miniconda3/etc/profile.d/conda.sh
conda activate phyluce-1.7.3

phyluce_assembly_match_contigs_to_probes \
    --contigs /export/home/a-e/amateos/phyluce_test_2/phylo_tests/ahya_tests/contigs \
    --probes /export/home/a-e/amateos/phyluce_test_2/data/hexa-v2-probes.fasta \
    --output uce_search_results_phylo_test