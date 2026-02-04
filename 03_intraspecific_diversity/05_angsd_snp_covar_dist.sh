#!/bin/bash

#12_angsd_snp_covar_dist_admix.sh located in 
#This script submits jobs to bin/angsd_ngs_dist.slurm with a set of parameters.
#Jobs are done with specified reference and sample set.

# Parameters for "angsd_ngs_dist.slurm" :
# $1=reffile $2=sampleinfo $3=bamdir $4=outdir
# (automated depths, minind; preset maxdepind and hwepval cutoffs.)


submit3(){
sbatch --partition=cpuq \
--job-name=snpAH20k \
--output=/export/home/a-e/amateos/phyluce_test_2/mapping_tests/logs/top100_uces_angsd_ngsdist_logs/angsd_$1.log \
  /export/home/a-e/amateos/phyluce_test_2/mapping_tests/angsd_ngs_dist.slurm $1 \
  /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/meta_file_$2.list \
  /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/10b-11_gatk_indel_angsd_depth/$2 \
  /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/top100_uces_12_ngs_dist_region/$2 $3 $4
}

#Full dataset run for 5-80%, maxDepthInd4, GL2, no-minHWEpval and uisng a region-filter file. ANGSD v0.929.


submit3 ahya_20k ahyaGRSamples 1e-3 1

# The "submit" function above calls the SLURM command "sbatch" for the job script
# "amgsd_ngs_dist.slurm" with the following parameters:

#'submit3' parameters
  #$1=reffile    (ahya500 ahya1k etc.)
  #$2=sample_set (ahyaSamples, etc.)
  #$3=snp_pval   (1e-3 etc.)
  #$4=prior      (1 or 2; -dopost)
