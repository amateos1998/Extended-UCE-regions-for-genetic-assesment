#!/bin/bash
#SBATCH --job-name=map2ref
#SBATCH --output=map2ref_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --partition=cpuq


#This script maps clean reads (PE and merged) to reference UCE references, then merges bam files into one after filtering steps.
#Prerequisite is to have bowtie2 databases already made under the names specified in 07-10a_map_pichard_samtools.sh'.
#This script is based on Therkildsen&Palumbi (2016) and Sato et al (2023). 

#when running manually with bash fro testing, set default values for slurm variables
####slrm_JOB_ID=${SLURM_JOB_ID:-map2ref}
#####slrm_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-4}

#log slurm info
echo "# Job ID: $SLURM_JOB_ID started"
echo "# CPUs available: $SLURM_CPUS_PER_TASK"
hostname
date
echo

#variables
sample=$1  
indir=$2
infolist=$3 #data/sample_info_test.list
refdir=$4   #data/ref
reflist=$5  #'ahya_500 ahya_1k etc'
outdir=$6   #results/07_map_reads
pdir=$(pwd)

echo "$sample reads are mapped"
echo "job ID: $SLURM_JOB_ID started"
hostname
date
echo

#Ensure /scratch is not full.
dir="/your/path/scratch"   
reqSpace=10     #in GByte
availSpace=$(df -BG "${dir}" | awk 'NR==2 { print $4+0 }')
if (( availSpace < reqSpace )); then
  echo "not enough Space" >&2
  exit 1
fi

#Copy files
echo "########## file copy to/your/path/mapping_tests/scratch"
mkdir -p /your/path/scratch/tmp.$SLURM_JOB_ID
cd /your/path/scratch/tmp.$SLURM_JOB_ID
cp /your/path/results/02-05_clean_reads/${sample}*.fastq.gz .
unpigz ${sample}*.fastq.gz
echo "# PE files and SE (merged) copied, unzipped"
ls -l *
cp $infolist info.list
echo

echo "extracting values from list"
#Extract relevant values from the sample list (sample ID, library ID, and platform unit)
 # SampleID	LibraryID	PlatformUnit	
SAMPLE_ID=$(grep -P "${sample}\t" info.list | cut -f 1)
LIB_ID=$(grep -P "${sample}\t" info.list | cut -f 2)
PU=$(grep -P "${sample}\t" info.list | cut -f 3)

echo "values extracted from list"
echo

#If multiple reference files are specificed, mapping steps will repeat for different files.
for reffile in $reflist; do   #Calling name of bt2 referernces
 cp $refdir/$reffile.*bt2 .   #copy bowtie2-library
 echo "############### ${reffile} #################"
 echo "########## bowtie2 started (paired)"
 bowtie2 -q --phred33 --very-sensitive-local -p $SLURM_CPUS_PER_TASK -I 0 -X 1500 --fr\
  --rg-id $SAMPLE_ID --rg SM:$SAMPLE_ID --rg LB:$LIB_ID --rg PU:$PU --rg PL:ILLUMINA -x $reffile\
  -1 ${sample}.notCombined_norm_1.fastq -2 ${sample}.notCombined_norm_2.fastq -S ${sample}_Paired_bowtie2_${reffile}.sam
 echo " Convert sam to bam file for storage"
 samtools view -bS -F 4 -@ $SLURM_CPUS_PER_TASK ${sample}_Paired_bowtie2_${reffile}.sam > ${sample}_Paired_bowtie2_${reffile}.bam
 rm ${sample}_Paired_bowtie2_${reffile}.sam
 echo
 echo "########## bowtie2 started (single)"
 bowtie2 -q --phred33 --very-sensitive-local -p $SLURM_CPUS_PER_TASK\
  --rg-id $sample --rg SM:$SAMPLE_ID --rg LB:$LIB_ID --rg PU:$PU --rg PL:ILLUMINA -x $reffile\
  -U ${sample}_SingleEndQualFiltered_norm.fastq -S ${sample}_SingleEnd_bowtie2_${reffile}.sam
 echo " Convert sam to bam file for storage"
 samtools view -bS -F 4 -@ $SLURM_CPUS_PER_TASK ${sample}_SingleEnd_bowtie2_${reffile}.sam > ${sample}_SingleEnd_bowtie2_${reffile}.bam
  # '-b to BAM-output, -S to auto input-format, -F 4 to exclude any unmapped (FLAGs coantain bit4) reads'
 rm ${sample}_SingleEnd_bowtie2_${reffile}.sam
 echo
 echo "########## bam files filtering and sorting"
 # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a mapping quality score < 20)
 samtools view -h -q 20 -@ $SLURM_CPUS_PER_TASK ${sample}_Paired_bowtie2_${reffile}.bam | grep -v XS:i\
  | samtools view -buS -@ $SLURM_CPUS_PER_TASK - | samtools sort -@ $SLURM_CPUS_PER_TASK - > ${sample}_paired_to_${reffile}_MinQ20_sorted.bam
 samtools view -h -q 20 -@ $SLURM_CPUS_PER_TASK ${sample}_SingleEnd_bowtie2_${reffile}.bam | grep -v XS:i\
  | samtools view -buS -@ $SLURM_CPUS_PER_TASK - | samtools sort -@ $SLURM_CPUS_PER_TASK - > ${sample}_SingleEnd_to_${reffile}_MinQ20_sorted.bam
 echo
 
 ## Filter the paired mappings by pairing status (separated these out to count the proportion of mapped pairs falling in each of these categories)
 # Extracting the concordantly mapped read pairs
 samtools view -h -b -f 0x002 -@ $SLURM_CPUS_PER_TASK ${sample}_paired_to_${reffile}_MinQ20_sorted.bam > ${sample}_paired_to_${reffile}_MinQ20_CP.bam
 # Extract reads where only one mate mapped (unpaired)
 samtools view -h -@ $SLURM_CPUS_PER_TASK ${sample}_paired_to_${reffile}_MinQ20_sorted.bam | perl -n -e 'print $_ if (/^\@/ || /YT:Z:UP/) '\
  | samtools view -bS -@ $SLURM_CPUS_PER_TASK - > ${sample}_paired_to_${reffile}_MinQ20_AllUP.bam
 # Could also extract unpaired reads through filtering on the bit flag (samtools view -bu -f 0x008),
 # but this filtering leaves out those pairs where one mate mapped with poor quality (and therefore got
 # filtered out so that the mapping was effectively unpaired)

 # Extract discordantly mapped read pairs
 samtools view -h -@ $SLURM_CPUS_PER_TASK ${sample}_paired_to_${reffile}_MinQ20_sorted.bam | perl -n -e 'print $_ if (/^\@/ || /YT:Z:DP/) '\
  | samtools view -bS -@ $SLURM_CPUS_PER_TASK - > ${sample}_paired_to_${reffile}_MinQ20_DiscordantlyMappedPairs.bam
 # Could also extract discordantly read pairs through filtering on the bit flag
 # (samtools view -b -F 0x0008 ${sample}_paired_to_${reffile}_MinQ20_sorted.bam | samtools view -b -F 0x002),
 # but that includes pairs where one mate mapped with poor quality (and therefore got filtered out so that the mapping was effectively unpaired)
 echo

 echo "########## Clipping overlaps"
 # Closer examination of the ${sample}_paired_To_NR_Trans_MinQ20_CP.bam revealed that there were a non-negligible number of
 # concordantly mapped read pairs for which the mapped ends overlapped (i.e. they had not been merged by flash).
 # To avoid double-counting these overlapping sections in SNP calling and genotype likelihood estimation, we used clipOverlap
 # to soft-clip overlapping read ends (maintaining only the read with the highest quality score in overlapping regions).
 
echo "######### check bam file is valid (no output=valid)#######"
  samtools quickcheck -v ${sample}_paired_to_${reffile}_MinQ20_CP.bam
echo

echo "###starting"
####forced to run on single thread, clipoverlap has issues in multithreaded envs. used TRUE flag to supress error
 ###taskset -c 0 bam clipOverlap --in ${sample}_paired_to_${reffile}_MinQ20_CP.bam \
  ####--out ${sample}_paired_to_${reffile}_MinQ20_CP_Clipped.bam \
  #####--storeOrig OC --stats --unmapped || true

####run clipoverlap######
bam clipOverlap --in ${sample}_paired_to_${reffile}_MinQ20_CP.bam\
  --out ${sample}_paired_to_${reffile}_MinQ20_CP_Clipped.bam\
  --storeOrig OC --stats --unmapped

# Stop script if output BAM file is not found, due to a fail in clipOverlap
if [ ! -f "${sample}_paired_to_${reffile}_MinQ20_CP_Clipped.bam" ]; then
    echo "Error: clipOverlap failed. Stopping script." >&2
    exit 1
fi

echo "### BAM clipping completed successfully."
 # Remove the reads that became unmapped in the clipping
 samtools view -hb -F 4 -@ $SLURM_CPUS_PER_TASK ${sample}_paired_to_${reffile}_MinQ20_CP_Clipped.bam\
  > ${sample}_paired_to_${reffile}_MinQ20_CP_Clipped_MappedOnly.bam
 echo

echo "########## bam files merging"
 # Merge the bam files to create a single bam file for each sample with all the reads to use for downstream analysis
 # (the overlap-clipped concordantly mapped reads, the unpaired mapped reads (where only one end of the pair mapped to
 # the transcriptome reference), and the single-end/orphaned mapped reads (the merged pairs and reads for which the other
 # mate in the pair had been lost during the quality filtering)). We did not include discordantly mapped pairs in any downstream analysis.
 samtools merge -f -@ $SLURM_CPUS_PER_TASK ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam\
  ${sample}_paired_to_${reffile}_MinQ20_CP_Clipped_MappedOnly.bam\
  ${sample}_paired_to_${reffile}_MinQ20_AllUP.bam\
  ${sample}_SingleEnd_to_${reffile}_MinQ20_sorted.bam
 echo "#####bam merged"
echo

echo "#####fixing BAM headres, picard does not accept the {} and ' characters"

samtools view -H ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam | grep '^@SQ' | awk '{print $2}' | head -n 5
samtools view -H ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam \
    | sed "s/[\{\}']//g" > fixed_header.sam
# Reheader BAM file
samtools reheader fixed_header.sam ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam > fixed.bam

# Replace the original file
mv fixed.bam ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam

samtools view -H ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam | grep '@SQ' | head -n 5

echo "#####file fixed"

 #Deduplication with picard
 echo "########## picard started"
/your/path/miniconda3/bin/picard MarkDuplicates I=${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam\
 O=${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup.bam M=${sample}_to_${reffile}_MinQ20_PairedAndSingle_dupstat.txt\
 VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
 echo "picard complete on $reffile"
 rm ${sample}_to_${reffile}_MinQ20_PairedAndSingle.bam
 echo

 echo "##### samtool sorting and indexing"
 samtools sort -o ${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam ${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup.bam -@ $SLURM_CPUS_PER_TASK
 samtools index ${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam -b
 echo "done"
 ls -l *
 echo

 echo "##### moving sorted bam and index files to $outdir/$reffile"
 mkdir -p $outdir/$reffile
 pigz ${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam*
 mv ${sample}_to_${reffile}_MinQ20_PairedAndSingle_dedup_sorted.bam*.gz $outdir/$reffile
 echo "done"
 date
 echo
 echo
done

rm -r /your/path/scratch/tmp.$SLURM_JOB_ID/*    #delete all remaining files in tmp folder.
echo "# job ID: $SLURM_JOB_ID finished."

