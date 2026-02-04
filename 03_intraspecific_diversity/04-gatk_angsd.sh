#!/bin/bash

#10b-11_indel_gatk_angsd_depth.sh
#This script submits a list of samples to bin/gatk_indel_angsd_depth.slurm
#Prerequisite is that all the mapped files are linked in to one folder specified below as $bamdir.


# The "submit" function below calls the SLURM command "sbatch" for the job script
# "bin/gatk_indel_angsd_depth.sbatch" with the following parameters:
#
#   $1 = reference file 
#   $2 = sample info
#   $3 = bamdir
#   $4 = outdir

submit(){
sbatch --partition=cpuq \
--job-name=gatkAHwg \
--output=/export/home/a-e/amateos/phyluce_test_2/mapping_tests/logs/ahya_GR_gatk/gatkAngsd_$3.log \
/export/home/a-e/amateos/phyluce_test_2/mapping_tests/gatk_indel_angsd_depth.slurm $1 \
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/data/meta_file_$2.list \
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/07_map_reads/$3 \
    /export/home/a-e/amateos/phyluce_test_2/mapping_tests/results/10b-11_gatk_indel_angsd_depth/$2 $4
}

#'submit' parameters
  #$1=reffile    (ahya_500 ahya_1k ahya_5k etc.)
  #$2=sample_set (meta_file_ahya)
  #$3=bamdir     (ahya_500 ahya_1k etc.)
  #$4=maxdepth   (~#sample x 10; eg. 1100)

submit ahya_WG ahyaGRSamples ahya_WG_run 1200