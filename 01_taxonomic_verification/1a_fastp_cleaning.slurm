#!/bin/bash
#SBATCH --job-name=
#SBATCH --output=
#SBATCH --partition=
#SBATCH --time=
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --array=1-8

#This script uses an array job in fastp to clean paired end reads

# Variables
READ_PREFIX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" read_names.txt) #txt file with read prefixes
INPUT_DIR="/your/path"
OUTPUT_DIR="/your/path"
JSON_OUT="/your/path/json_html"
HTML_OUT="/your/path/json_html"

R1="${INPUT_DIR}/${READ_PREFIX}_1.fq.gz"
R2="${INPUT_DIR}/${READ_PREFIX}_2.fq.gz"
OUT_R1="${OUTPUT_DIR}/${READ_PREFIX}.out.R1.gz"
OUT_R2="${OUTPUT_DIR}/${READ_PREFIX}.out.R2.gz"

# Debugging logs
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "READ_PREFIX=${READ_PREFIX}"
echo "Input file R1: $R1"
echo "Input file R2: $R2"
echo "Output file OUT_R1: $OUT_R1"
echo "Output file OUT_R2: $OUT_R2"

# Check input files
if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Input files missing: $R1 or $R2"
    exit 1
fi

# Run fastp
fastp -i "$R1" -I "$R2" -o "$OUT_R1" -O "$OUT_R2" \
    --json "${JSON_OUT}/${READ_PREFIX}.json" \
    --html "${HTML_OUT}/${READ_PREFIX}.html" 2>&1 | tee -a /your/path
