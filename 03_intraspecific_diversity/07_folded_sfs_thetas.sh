#!/bin/bash

#18_folded_sfs_thetas.sh located in /opt/extern/bremen/symbiosis/sato/05_exome_SNPs_Oalg_180115/
#This script submits a list of samples from a populaiton to bin/sfs_folded_thetas.qsub
#Samples of each groups are listed in data/sample_info_180312_SOR.list etc.

#sample_info_XXX.list contains, without headers and tab-delimiated:
#SampleID     LibraryID       PlatformUnit    Haplotype       Location        ReadLength
#OalgA10CA 2162_G_1        MPIPZ           A               Cavoli_Elba     150
#OalgA10SA 2162_Q_1        MPIPZ           A               SantAndrea_Elba 150
#OalgA1CA  CSP2012_1049310 JGI             A               Cavoli_Elba     150
#OalgA1SA  CSP2012_1021950 JGI             A               SantAndrea_Elba 150

# ALL MAPPED RESULTS WERE LINKED TO ONE FOLDRER $BAMDIR!!*********


###################################################################################################
submit_folded_mI50(){
sbatch --job-name=thetasAh20k \
--output=/export/home/a-e/amateos/phyluce_test_2/mapping_tests/logs/top100_uces_18thetassfs/theta_$1.log \
 /export/home/a-e/amateos/phyluce_test_2/mapping_tests/sfs_folded-only_thetas_minInd50.slurm $1 $2\
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/top_uces_18_folded_sfs_thetas_minInd50\
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/10b-11_gatk_indel_angsd_depth/$3\
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/ahya_gr_meta.list

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






