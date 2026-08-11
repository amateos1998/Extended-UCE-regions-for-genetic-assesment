# Extended-UCE-regions-for-genetic-assesment
### Input files required to run the pipline
You will need genomic paired end reads in fasta format (fq.gz), a UCE baitset designed for your study taxa (in this case we used the hexa-v2 baitset develop by [Cowman et al. 2020](https://www.sciencedirect.com/science/article/pii/S1055790320302165/)), and alignments from an existing phylogenetic tree to place your samples into (fasta or nexus). 

You will also need reference genomes to generate the extended UCE references. You can generate a UCE reference using a de-novo assembly of one of your lcWGS short reads, however this will limit the maximum size of the UCE flanking region to the length of your contigs (we recommend using published genomes if possible).

It's also a good idea to have a master metadata file (.list or .txt) with all the samples and to have metadata files for each group of study (by population, species, haplotype, etc.)

Note: lcWGS ilumina paired-end short reads are best in this pipeline, depth and coverage of 4-5x is recommended for optimal results.

### Software requirements to run the pipeline
The whole pipeline is based in Linux using bash and R.
You will need to have access to a HPC workload manager (in this case I used SLURM, but you can use qsub or any other manager).

You will need to install the following tools to run the pipeline (I recomend using conda to manage all your tools, all my scripts use conda to activate different environments):
* Phyluce v1.7.3 (the phyluce environment already comes with some of the these binaries)
* fastp v0.23.2
* MAFFT v7.0
* IQ-TREE v2.3.6
* SPAdes Genome Assembler v4.1.0
* fastuniq v1.1
* Trimmomatic v.0.33
* bbsplit v.35.85
* FLASH v1.2.11
* bowtie2 v2.5.4
* faToTwoBit
* GATK v3.8b
* ANGSD v0.928
* ngsDist v1.0.1
* cdhit v4.8.1

You will also need some standard tools like gzip to work with all the files. It's important to use the versions specified in this list, as some features used in this pipeline might not be available in newer builds.

Step 03_intraspecific_diversity is bases on X
