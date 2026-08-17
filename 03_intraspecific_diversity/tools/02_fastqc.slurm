#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=24:00:00
#SBATCH --partition=cpuq

#This script performs FastQC



#log slurm info
echo "# Job ID: $SLURM_JOB_ID started"
echo "# CPUs available: $SLURM_CPUS_PER_TASK"
hostname
date
echo

# Variables (three input files and output directory)
infile1=$1
infile2=$2
infile3=$3
outdir=$4

echo "# Job ID: $SLURM_JOB_ID started"
hostname
date
echo "Processing files:"
echo "$infile1"
echo "$infile2"
echo "$infile3"
echo

# Ensure /scratch is not full
dir="/your/path/scratch"   
reqSpace=10     # in GByte
availSpace=$(df -BG "${dir}" | awk 'NR==2 { print $4+0 }')
if (( availSpace < reqSpace )); then
  echo "Not enough space on scratch!" >&2
  exit 1
fi

# Copy only the input files to scratch
echo "########## Copying files to /your/path/scratch"
mkdir -p /your/path/scratch/tmp.$SLURM_JOB_ID
cd /your/path/mapping_tests/scratch/tmp.$SLURM_JOB_ID

cp "$infile1" .
cp "$infile2" .
cp "$infile3" .
sync  # Ensure copy is complete before proceeding

echo "# Files copied:"
ls -l "$(basename "$infile1")"
ls -l "$(basename "$infile2")"
ls -l "$(basename "$infile3")"
echo

# Run FastQC on all three files
echo "Running FastQC on all three files..."
fastqc --threads "$SLURM_CPUS_PER_TASK" --nogroup -o "$outdir" "$(basename "$infile1")" "$(basename "$infile2")" "$(basename "$infile3")"

echo "FastQC completed successfully."

# Move results back (restored original commented-out lines)
#echo "moving result"
#mkdir -p $outdir
#mv *fastqc.zip *fastqc.html $outdir

# Clean up scratch directory
rm -r /your/path/scratch/tmp.$SLURM_JOB_ID/*

echo "Job ID: $SLURM_JOB_ID finished."
date
