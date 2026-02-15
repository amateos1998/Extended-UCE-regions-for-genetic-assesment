#!/bin/bash
#SBATCH --job-name=gatk
#SBATCH --output=gatkAhya_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=9-48:00:00
#SBATCH --partition=cpuq

#gatk_indel_angsd_depth.qsub located in bin
#This script takes sorted bam files and reference file and realigns around intron/exon boundaries for better SNP detection,
# followed by depth-identification by ANGSD for fiding out depth parameters in the next step.
#Samples included are listed in $2 sampleinfo.

#log slurm info
echo "# Job ID: $SLURM_JOB_ID started"
echo "# CPUs available: $SLURM_CPUS_PER_TASK"
hostname
date
echo

reffile=$1 #Reference file
sampleinfo=$2 #data/x_sample.list
bamdir=$3 #results/07-10a_mapping
outdir=$4 #results/10b-11_indel_angsd_depth
pdir=$(pwd)
maxdepth=$5 #1100 etc. (~#sample x 10. 100 samples per spp.)

echo "job ID: $SLURM_JOB_ID started"
hostname
date
echo

#Ensure /scratch is not full.
dir="/export/home/a-e/amateos/phyluce_test_2/mapping_tests/scratch"   
reqSpace=50     #in GByte
availSpace=$(df -BG "${dir}" | awk 'NR==2 { print $4+0 }')
if (( availSpace < reqSpace )); then
  echo "not enough Space" >&2
  exit 1
fi

#Copy files

mkdir -p /export/home/a-e/amateos/phyluce_test_2/mapping_tests/scratch/tmp.$SLURM_JOB_ID
cd /export/home/a-e/amateos/phyluce_test_2/mapping_tests/scratch/tmp.$SLURM_JOB_ID




##### REFERENCE FILE HERE #####

echo "####reference files copying#####"
cp /export/home/a-e/amateos/phyluce_test_2/mapping_tests/reference/ahya_WG.fasta $reffile.fasta





echo
echo "###reference copied#####"
echo

echo "############### ${reffile} #################"
echo "##### copying all sorted bam and index files, unzipping them, and making a list of samples included"
echo "Sample list: $sampleinfo"
echo "Source: $bamdir/$reffile"
for sample in $(cat $sampleinfo | cut -f1); do
 echo $sample
 cp $bamdir/$reffile/${sample}_*bam.gz .
done
unpigz *bam.gz
echo

for sample in $(cat $sampleinfo | cut -f1); do
 cp $bamdir/$reffile/${sample}_*bai.gz .
done
unpigz *bai.gz

######removed the | grep "Oalg" flag

ls *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam > Samples_To_${reffile}_dedup_sorted_bam.list
echo "done"
echo
echo
echo
echo "#####picard doesn't like {} and ', fixing file headers for next step#####"
echo
echo "##### Fixing reference file headers (removing {} and other invalid characters)"
sed -i "s/[\{\}']//g" ${reffile}.fasta

echo "##### Verifying sequence names"
grep ">" ${reffile}.fasta | grep -E "[^0-9A-Za-z!#$%&+./:;?@^_|~-]" | head -n 5 && echo "Warning: Some invalid characters remain!" || echo "Headers look good."


echo "##### Process completed successfully."
echo "##### File fixed."
echo
echo
echo

echo "##### indexing reference file and make dictionary file"
samtools faidx $reffile.fasta
/export/home/a-e/amateos/miniconda3/bin/picard CreateSequenceDictionary \
 R=${reffile}.fasta O=${reffile}.dict
echo "done"
echo

echo "########SWITCHING TO JAVA 1.8 FOR GATK 3.8##########"
# Save current JAVA_HOME and PATH
ORIGINAL_JAVA_HOME=$JAVA_HOME
ORIGINAL_PATH=$PATH

# Set Java 1.8 (for GATK)
export JAVA_HOME=/export/home/a-e/amateos/tools/java/jdk8u392-b08
export PATH=$JAVA_HOME/bin:$PATH
# Confirm Java version (should be Java 1.8)
java -version
echo
echo
echo "##### gatk identifiying target intervals"
java -Xmx250G -jar /export/home/a-e/amateos/tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK3.jar \
  -T RealignerTargetCreator \
  -R $reffile.fasta \
  -I Samples_To_${reffile}_dedup_sorted_bam.list \
  -o ${reffile}_indel.intervals \
  -drf BadMate \
  -nt $SLURM_CPUS_PER_TASK
echo "done"
echo

echo "##### gatk re-aligning around indels"
java -Xmx200G -jar /export/home/a-e/amateos/tools/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK3.jar \
   -T IndelRealigner \
   -R $reffile.fasta \
   -I Samples_To_${reffile}_dedup_sorted_bam.list \
   -targetIntervals ${reffile}_indel.intervals \
   --consensusDeterminationModel USE_READS \
   --nWayOut _realigned.bam
echo "done"
echo
# Restore original Java (Miniconda Java)
export JAVA_HOME=$ORIGINAL_JAVA_HOME
export PATH=/export/home/a-e/amateos/miniconda3/bin:$ORIGINAL_PATH
echo "#######restoring java to version before GATK######"
# Confirm restoration (should return Miniconda Java path)
java -version
echo
echo "##### copying results of indel-align"
echo "Target: $outdir/$reffile"
mkdir -p $outdir/$reffile
rm *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam.bai
for bamfile in *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted_realigned.bam; do pigz -c $bamfile > $bamfile.gz; done
for baifile in *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted_realigned.bai; do pigz -c $baifile > $baifile.gz; done
mv *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted_realigned.bam.gz\
 *_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted_realigned.bai.gz $outdir/$reffile
echo

echo "##### starting to calculate reference mapping depth with ANGSD"
ls *realigned.bam > bam_list.txt
n_ind=$(cat bam_list.txt | wc -l)          #Calculate the number of individual samples
echo "BAM list containing $n_ind samples is made"
echo
echo "ANGSD starting to check the mapping depth profile..."
date
samtools faidx $reffile.fasta
/export/home/a-e/amateos/phyluce_test_2/tools/ngsTools/angsd/angsd -P $SLURM_CPUS_PER_TASK -b bam_list.txt -out ${n_ind}Samples_to_${reffile}_MinQ20_qc -ref $reffile.fasta\
 -minQ 2 -trim 0 -remove_bads 0 -uniqueOnly 0 -only_proper_pairs 0 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth $maxdepth
Rscript /export/home/a-e/amateos/phyluce_test_2/tools/ngsTools/Scripts/plotQC.R ${n_ind}Samples_to_${reffile}_MinQ20_qc 2> /dev/null
 # Use -only_proper_pairs 0 to also retain reads for which only one end of a pair is mapped (because mapping to a transcriptome, not the whole genome)
 # Remove_bads 0 and -uniqueOnly 0 as they were already filtered from the bam files for poorly mapped reads.

echo "ANGSD complete"
echo date
echo

echo "########## moving result to:"
echo " $outdir/$reffile"
rm *_MinQ20_PairedAndSingle_dedup_sorted_realigned.ba? bam_list.txt *fasta* *bam.list *qc.qs *intervals *dict   #(optional) remove intermediate files not used later.
mv * $outdir/$reffile
rm -r /scratch/tmp.$JOB_ID
echo
echo

echo "job ID: $JOB_ID finished."
date
