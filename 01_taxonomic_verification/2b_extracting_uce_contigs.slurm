#!/bin/bash
#SBATCH --job-name=extracting_uce_contigs     # Job name
#SBATCH --partition=cpuq   # Specify the partition you want to use
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --time=48:00:00                   # Time limit hrs:min:sec
#SBATCH --mem=20G                          # Memory per node
#SBATCH --output=/your/path/logs/extractUCE_%j.out  # Standard output log

###THIS SCRIPT extracts UCE loci from contigs assembled from lcWGS

# Activate the conda environment where Phyluce is installed
source /your/path/miniconda3/etc/profile.d/conda.sh
conda activate phyluce-1.7.3


#This set of code takes the .lastz files from finding UCEs and turns them into .fasta files for each locus and sample in
#the taxon set (ie. config file)

date
echo "get fastas from match counts"


phyluce_assembly_get_fastas_from_match_counts \
    --contigs /your/path \
    --locus-db /your/path/uce_search_results/probe.matches.sqlite \
    --match-count-output /your/path/uces_extracted_contigs/all-taxa-incomplete.conf \
    --output //your/path/all-taxa-incomplete.fasta \
    --incomplete-matrix /your/path/uces_extracted_contigs/all-taxa-incomplete.conf



echo "explode the monolithic FASTA by locus (can also do by --taxa)"
phyluce_assembly_explode_get_fastas_file \
    --input all-taxa-incomplete.fasta \
    --output exploded-fastas 

# get summary stats on the FASTAS
for i in exploded-fastas/*.fasta;
do
    phyluce_assembly_get_fasta_lengths --input $i --csv;
done

echo "extraction step 2 done"


#optional step: you can add other incomplete fastas from different runs of step 2b at this point to add new samples to your tree. Simply concatenate all desired fasta files into one.

#phyluce_align_concatenate_alignments --alignments top100_uces_aligs --input-format #fasta --output top100_uce_tree --phylip
	
