# Extended-UCE-regions-for-genetic-assesment

### Requirements to run the pipeline
The whole pipeline is based in Linux using bash and R.
You will need to have access to a HPC workload manager (in this case I used SLURM, but you can use qsub or any other manager)
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
* GATK v3.8
* ANGSD v0.928
* ngsDist v1.0.1
* cdhit v4.8.1

You will also need some standard tools like gzip to work with all the files
