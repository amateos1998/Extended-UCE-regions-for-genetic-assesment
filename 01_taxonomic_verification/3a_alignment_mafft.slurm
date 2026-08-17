#!/bin/bash

#SBATCH --job-name=alUCEs
#SBATCH --time=99:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=250G
#SBATCH --partition=cpuq
#SBATCH --output=/your/path/



#This script aligns your uce alignments to the reference uce alignments (Cowman et al. 2020) locus by locus.

# Load MAFFT through Phyluce
source /your/path/miniconda3/etc/profile.d/conda.sh
conda activate phyluce-1.7.3

date
echo

for i in *.unaligned.fasta; do
outname=$(basename "$i" .unaligned.fasta)

echo "Aligning $i to $outname.fasta"

mafft --auto --adjustdirection \
 --add "$i" \
 --reorder --thread 30 "$outname.fasta" \
 > /your/path/uces_extracted_contigs/uce-loci/${outname}.fasta;

echo "Aligned $i to ${outname}.fasta, saved to uce-loci/${outname}.fasta"
echo

done

date
echo

