#!/bin/bash
#SBATCH --job-name=top75pTree
#SBATCH --output=/export/home/a-e/amateos/phyluce_test_2/final_tree/results/uces_extracted_contigs/logs/top100_loci_75p_%j.log
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=20-24:00:00
#SBATCH --partition=cpuq


date
echo
echo "starting IQtree"
echo

/export/home/a-e/amateos/phyluce_test_2/tools/iqtree-2.3.6-Linux-intel/bin/iqtree2 -p top100-trimmed-75p \
	--prefix topconcat75edge \
	-B 1000 \
	--boot-trees \
	-T auto \
	--threads-max 10 \
	-m TESTMERGE \
	-alninfo \
	--merge-model GTR \
	--merge-rate G \
	--rcluster 10 \
	--cptime 4000


echo
echo "iqtree done"
echo 
date