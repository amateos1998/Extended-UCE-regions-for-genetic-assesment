#!/bin/bash
#SBATCH --job-name= #Job name
#SBATCH --partition= #Specify the partition you want to use (cpuq, memq, etc)
#SBATCH --cpus-per-task= #Number of CPU cores per task
#SBATCH --time= #Time limit hrs:min:sec
#SBATCH --mem= #Memory per node
#SBATCH --output= #Standard output log


###This script harvests UCE loci from genomes, before running, convert fastas to .2bit ###(you can use the faToTwoBit tool). Create a folder named after the genome file with ###the .tab and .2bit files.

#The script uses a genomes.conf file which should have the following format:
#[scaffolds]
#species_a:/your/path/species_a.2bit
#species_b:/your/path/species_b.2bit
#...
#species_z:/your/path/speciesz.2bit

date
echo

# Activate the conda environment where Phyluce is installed
source /your/path
conda activate phyluce-1.7.3

echo "phyluce searching for UCEs"

phyluce_probe_run_multiple_lastzs_sqlite \
    --db \
    --output \
    --scaffoldlist \ #list of genomes
    --genome-base-path  \
    --probefile /your/path/Cowman_etal_APPENDIX_D-hexa-v2-sclerac-subset-final-probes.fasta \
    --cores 12

date 
echo "UCE search done"
