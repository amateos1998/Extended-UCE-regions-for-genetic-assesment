#!/bin/bash
#SBATCH --job-name=ahya_asmbly          # Job name
#SBATCH --partition=cpuq                    # Partition name (adjust according to your cluster)
#SBATCH --time=99:00:00                     # Time limit (hh:mm:ss)
#SBATCH --cpus-per-task=60                  # Number of CPUs per task
#SBATCH --mem=500G                          # Memory allocation
#SBATCH --output=/export/home/a-e/amateos/phyluce_test_2/logs/spades_assembly_%j.log     # Standard output log (%j = job ID)

#This script performs de-novo assembly of your samples for downstream analysis in phyluce
date
echo

# Activate Conda and set up the environment
source /your/path
conda activate spades_env 


# Directories
READS_DIR="/your/path"
OUT_DIR="/your/path"

# Number of cores
THREADS=60

# List of sample IDs
SAMPLES=("sample_x")



# Loop over each sample and run SPAdes for each one
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"

    # Input FASTQ files for the sample
    R1="$READS_DIR/${SAMPLE}_R1.fq.gz"
    R2="$READS_DIR/${SAMPLE}_R2.fq.gz"
    
    # Output directory for the current sample
    SAMPLE_OUT_DIR="$OUT_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT_DIR"

    # Run SPAdes assembly
    spades.py \
        -1 "$R1" \
        -2 "$R2" \
        -o "$SAMPLE_OUT_DIR" \
        --threads $THREADS \
        --memory 500

    echo "Assembly for $SAMPLE completed."
done

echo "All assemblies completed."

echo "assembly QC"

# Define the output CSV file
output_csv="/your/path/assemblyQC.csv"

# Loop through all .fasta files in the specified directory and run phyluce_assembly_get_fasta_lengths
for i in your/path/*.fasta; 
do
    phyluce_assembly_get_fasta_lengths --input "$i" --csv >> $output_csv
done

echo "All FASTA file lengths written to $output_csv"

