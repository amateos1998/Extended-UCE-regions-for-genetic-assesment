#!/bin/bash

#18_folded_sfs_thetas.sh 
#This script submits a list of samples from a populaiton to bin/sfs_folded_thetas.qsub
#Samples of each groups are listed in data/sample_info_180312_SOR.list etc.




###################################################################################################
submit_folded_mI50(){
sbatch --job-name=thetasAh20k \
--output=/your/path/s/theta_$1.log \
 /your/path/sfs_folded-only_thetas_minInd50.slurm $1 $2\
   /your/path/results/top_uces_18_folded_sfs_thetas_minInd50\
   /your/path/results/10b-11_gatk_indel_angsd_depth/$3\
   /your/path/data/ahya_gr_meta.list

#Parameters for "bin/sfs_folded_thetas.qsub":
 #reffile=$1 #ahya_500, ahya_1k, etc.
 #pop=$2 #TAKEN FROM META FILE
 #outdir=$3 #/results/18_folded_sfs_thetas
 #bamdir=$4 #/results/10b-11_indel_angsd_depth/
 #sampleinfo=$5 #/data/sample_info_ahya
}

submit_folded_mI50 ahya_20k CBHE ahyaGRSamples



#'submit' parameters
  #$1=reffile 
  #$2=pop    
  #$3=bamdir  






