#!/bin/bash
#SBATCH --job-name=sfs_thetas
#SBATCH --output=sfs_thetas_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=9-48:00:00
#SBATCH --partition=cpuq



#This script makes ANGSD-saf files (.gz & .idx) and calculates folded SFS and thetas/D-stats using realSFS (unfolded SFS survey is skipped).
#SAF file is folded for thetaStat and generatd with ANGSD -doSaf 1 -fold 1. (not using -doSaf 2 inbreed version as coefficient was all 0).
#ANGSD v0.929 was used instead of v0.921 because of a fatal error fixed for thetas. thetaStat and realSFS are also from ANGSD 0.929.
#Uses a minInd 50% filter and region filter file (see below).

pdir=$(pwd)
reffile=$1 #Oalg_verC3 etc.
pop=$2 #MAGL_R20 etc.
outdir=$3 #/results/18_folded_sfs_thetas
bamdir=$4 #/results/10b-11_gatk_indel_angsd_depth/aspatSamples/aspat_2k/ ETC..
sampleinfo=$5 #/data/meta_file_aspatSamples.list etc.

minind=$(cat $sampleinfo | wc -l | awk '{print int($1/2+0.5)}')   #10 etc. 50% of sample number, rounded up integral
upperdepth=$(cat $sampleinfo | wc -l | awk '{print $1*10}') #200 etc. 10X of #samples x (read-depths; e.g. 1X) as the upper limit of global depth-graph. Used for a depth histogram.
maxdepth=$(cat $sampleinfo | wc -l | awk '{print $1*20}') #200 etc. 10X of sample number for 1X-read depth data. Increase according to the read depth.
                                                          #Upper limit of depth range for ANGSD to exclude reads from repeat regions but accommodate stochasticity.

### REGION FILE HERE ####

regions=/your/path/a_rf.txt

#log slurm info
echo "# Job ID: $SLURM_JOB_ID started"
echo "# CPUs available: $SLURM_CPUS_PER_TASK"
hostname
date
echo

#Ensure /scratch is not full.
dir="/your/path/scratch"   
reqSpace=50     #in GByte
availSpace=$(df -BG "${dir}" | awk 'NR==2 { print $4+0 }')
if (( availSpace < reqSpace )); then
  echo "not enough Space" >&2
  exit 1
fi

mkdir -p /your/path/scratch/tmp.$SLURM_JOB_ID
cd /your/path/scratch/tmp.$SLURM_JOB_ID

echo "copying reference file $reffile"
cp /your/path/reference/ahya_20k.fasta $reffile.fasta


echo "copying realingned bam and index files, unzipping them, and making a list of samples included"
echo "Population: $pop"; mkdir $pop
echo "Source: $bamdir/$reffile"
for sample in $(cat $sampleinfo | cut -f1); do
 echo $sample
 cp $bamdir/$reffile/$sample*bam.gz $pop/.
done
unpigz $pop/*bam.gz
for sample in $(cat $sampleinfo | cut -f1); do
 cp $bamdir/$reffile/$sample*bai.gz $pop/.
done
unpigz $pop/*bai.gz
ls $pop/*realigned.bam > bam_list.txt
n_ind=$(cat bam_list.txt | wc -l)          #Count the number of individual samples
echo
echo "done"
echo " "
#ls -l *
echo

echo "folded SFS calculation for thetas stats..."
echo "Calculating SAF (sample allele frequency) for $pop for folded SFS"
samtools faidx $reffile.fasta
/your/path/tools/angsd928/angsd -P $SLURM_CPUS_PER_TASK -b bam_list.txt -ref $reffile.fasta -anc $reffile.fasta -rf $regions \
 -out ${pop}_${n_ind}samples_to_${reffile}_folded -minInd $minind\
 -GL 2 -doCounts 1 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 0\
 -minQ 2 -setMaxDepth $maxdepth -doSaf 1 -fold 1 -underFlowProtect 1
  #-fold 1 - used to make folded SFS as the ancestory state is unknown for nucleotide diversity and neutuarity tests.
  #-underFlowProtect 1: For large data sets (large number of individuals) underflow projection is needed.
  #-GL 2; GATK genotype likelihood model is used (instead of 1; SAMtools) here.
echo "SAF for $pop complete"
date
echo
ls -l ${pop}_${n_ind}samples_to_${reffile}_folded.*
echo

echo "SFS estimation"
/your/path/tools/angsd928/misc/realSFS ${pop}_${n_ind}samples_to_${reffile}_folded.saf.idx -P $SLURM_CPUS_PER_TASK 2> /dev/null > ${pop}_${n_ind}samples_to_${reffile}_folded.sfs
echo "done"
ls -l ${pop}_${n_ind}samples_to_${reffile}_folded.sfs
cat ${pop}_${n_ind}samples_to_${reffile}_folded.sfs
echo



#using newer angs v.
echo "Calculating thetas for $pop using folded SFS"
samtools faidx $reffile.fasta
/your/path/tools/angsd928/angsd -P $SLURM_CPUS_PER_TASK -b bam_list.txt -ref $reffile.fasta -anc $reffile.fasta -rf $regions \
 -out ${pop}_${n_ind}samples_to_${reffile}_folded_nucdiv -minInd $minind\
 -GL 2 -doCounts 1 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 0\
 -minQ 2 -setMaxDepth $maxdepth -doSaf 1 -fold 1 -underFlowProtect 1\
 -doThetas 1 -pest ${pop}_${n_ind}samples_to_${reffile}_folded.sfs
  
#-fold 1 - used to make folded SFS, as the ancestory state is unknown for nucleotide diversity and neutuarity.
  #-underFlowProtect 1: For large data sets (large number of individuals) underflow projection is needed.
  #-doThetas 1 : to calculate statistics, -pest file.sfs to set the SFS as prior information.
  #angsd v0.929 was used instead of 0.921 due to a segmentation fault error.
  #-GL 2; GATK genotype likelihood model is used (instead of 1; SAMtools) here.

echo "Calculation thetas for $pop complete"
date
rm -rf $pop
echo
ls -l *
echo

echo "Nucleotide diversity estimation and neutuarity tests"


/your/path/tools/angsd928/misc/thetaStat do_stat ${pop}_${n_ind}samples_to_${reffile}_folded_nucdiv.thetas.idx \
 -outnames ${pop}_${n_ind}samples_to_${reffile}_folded_nucdiv.thetas #results are per chromosome.



echo

echo "all processes done."; echo
date
rm *.gz *.saf.idx bam_list.txt $reffile.fasta $reffile.fasta.fai
ls -l *
echo
echo "moving results to:"
echo " $outdir/$reffile"
mkdir -p $outdir/$reffile
mv * $outdir/$reffile
echo " "
rm -r /scratch/$USER/tmp.$SLURM_JOB_ID

echo "job ID: $SLURM_JOB_ID
 finished."
date
sleep 30
mv $pdir/$SLURM_JOB_NAME.o$SLURM_JOB_ID $outdir/$reffile/${pop}_to_${reffile}_sfs_folded_thetas.qsub_log
