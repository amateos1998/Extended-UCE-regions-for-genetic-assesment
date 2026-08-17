#!/bin/bash
#SBATCH --job-name=seqprep
#SBATCH --output=seqprep_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=9-24:00:00
#SBATCH --partition=cpuq

#bbsplit_norm_clean.qsub located in bin/
#This script takes raw paired fq.gz files and does de-duplication (fastuniq), adapter-clipping (trimmomatic)
#filtering out Oalg symbionts and mito DNA (bbsplit), QC-triming (trimmomatic), pair-read merging (flash), and decontamination (bowtie2).
#Scripts were partially based on Therkildsen&Palumbi (2016)


#log slurm info
echo "# Job ID: $SLURM_JOB_ID started"
echo "# CPUs available: $SLURM_CPUS_PER_TASK"
hostname
date
echo

#variables

sample=$1
outdir=$2
target_size=2000000000 # = 2.0 Gbases etc. Upper file size cap to avoid biassed noise. 
pdir=$(pwd)

echo "# job ID: $SLURM_JOB_ID started"
hostname
date
echo
echo "$sample processed."
echo "Sample-size target is $(awk -v tar=$target_size "BEGIN {printf \"%.1f\n\", tar / 1000000000}") Gbases of host nuc-DNA reads."

#Ensure /scratch is not full.
dir="/export/home/a-e/amateos/phyluce_test_2/mapping_tests/scratch"   
reqSpace=10     #in GByte
availSpace=$(df -BG "${dir}" | awk 'NR==2 { print $4+0 }')
if (( availSpace < reqSpace )); then
  echo "not enough Space" >&2
  exit 1
fi

#Copy files
echo "########## file copy to /your/path/scratch"
mkdir -p /your/path/scratch/tmp.$SLURM_JOB_ID
cd /your/path/scratch/tmp.$SLURM_JOB_ID
cp /your/path/data/Ahya_samples/${sample}*.fq.gz .
unpigz ${sample}*.fq.gz
echo "# PE files copied, unzipped"
ls * -l
echo

#Deduplication
echo "########## fastuniq-deduplicating..."
ls ${sample}*.fq > seqfile.list     #make an input-file list for fastuniq
fastuniq -i seqfile.list -t q -o ${sample}_FastUniq_1.fastq -p ${sample}_FastUniq_2.fastq -c 0
rm -f  seqfile.list ${sample}*.fq
echo "# fastuniq complete"
ls -l *
echo

echo "Job ID: '$SLURM_JOB_ID'"

echo "########## Adapter-clipping by trimmomatic..."
#Adapter-clipping
/your/path/miniconda3/envs/phyluce-1.7.3/bin/trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
 /your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_FastUniq_1.fastq \
 /your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_FastUniq_2.fastq \
/your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_AdapterClipped_forward_paired.fastq \
/your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_AdapterClipped_forward_unpaired.fastq \
/your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_AdapterClipped_reverse_paired.fastq \
 /your/path/scratch/tmp.${SLURM_JOB_ID}/${sample}_AdapterClipped_reverse_unpaired.fastq \
 ILLUMINACLIP:/your/path/bbmap/resources/adapters.fa:2:30:10:4:true/adapters.fa:2:30:10:4:true \
 -phred33
 # This specifies; seed mitmatch=2; palindrome clip threshold=30; simple clip threshold=10; min adapter length=4.

rm -f ${sample}_AdapterClipped_forward_unpaired.fastq ${sample}_AdapterClipped_reverse_unpaired.fastq
echo "# adapter clipping complete"
ls -l *
echo

#BBsplit into Oalg symbionts, mitochondria, and the rest containing host nucDNA reads, Omitted mito.
echo "########## BBsplit to Symb, Mito, and other(clean)"
echo "# ref file copying"
cp /your/path/data/zoox_genomes/Zoox_combined.fasta ./symbiont_genomes.fasta ##zoox genome
gunzip -c /your/path/bbmap/resources/phix174_ill.ref.fa.gz > phix.fasta

echo "bbsplit processing..."
bbsplit.sh in1=${sample}_AdapterClipped_forward_paired.fastq in2=${sample}_AdapterClipped_reverse_paired.fastq\
 ref=symbiont_genomes.fasta,phix.fasta basename=bbsplit_%_${sample}_#.fq\
 outu1=paired_clean1.fastq outu2=paired_clean2.fastq threads=$SLURM_CPUS_PER_TASK nodisk
rm -f ${sample}_AdapterClipped_forward_paired.fastq ${sample}_AdapterClipped_reverse_paired.fastq
echo
ls -l bbsplit_*fq paired_clean?.fastq
rm bbsplit_phix*.fq
pigz  bbsplit_*fq
echo


echo "########## QC-trimming by trimmomatic..."
#QC-trimming
/your/path/miniconda3/envs/phyluce-1.7.3/bin/trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 \
 /your/path/scratch/tmp.${SLURM_JOB_ID}/paired_clean1.fastq \
/your/path/scratch/tmp.${SLURM_JOB_ID}/paired_clean2.fastq \
 ${sample}_AdapterClipped_QualFiltered_forward_paired.fastq \
 ${sample}_AdapterClipped_QualFiltered_forward_unpaired.fastq \
 ${sample}_AdapterClipped_QualFiltered_reverse_paired.fastq \
 ${sample}_AdapterClipped_QualFiltered_reverse_unpaired.fastq \
 LEADING:10 SLIDINGWINDOW:4:2 MINLEN:50 #first 10-bases off if below q-threshold, window of 4 bases to check quality >2, shorter than 50 to dump.
rm -f paired_clean1.fastq paired_clean2.fastq
echo "# trimmomatic complete"
ls -l *
echo

echo "########## read-merging (flash) started"
flash -M 110 ${sample}_AdapterClipped_QualFiltered_forward_paired.fastq ${sample}_AdapterClipped_QualFiltered_reverse_paired.fastq\
 --output-prefix $sample -t $SLURM_CPUS_PER_TASK
rm -f ${sample}_AdapterClipped_QualFiltered_forward_paired.fastq ${sample}_AdapterClipped_QualFiltered_reverse_paired.fastq
echo "# flash complete"
ls -l *
echo

#Merge orphaned reads from trimming and flash-combined (extended) reads together as all single-end reads.
echo "rescuing QC-orphaned reads into flash-merged reads..."
cat ${sample}.extendedFrags.fastq\
 ${sample}_AdapterClipped_QualFiltered_forward_unpaired.fastq\
 ${sample}_AdapterClipped_QualFiltered_reverse_unpaired.fastq > ${sample}_SingleEndQualFiltered.fastq
rm -f ${sample}.extendedFrags.fastq\
 ${sample}_AdapterClipped_QualFiltered_forward_unpaired.fastq\
 ${sample}_AdapterClipped_QualFiltered_reverse_unpaired.fastq
ls -l *
echo


echo "#SKIP DECONTAMINATION STEP AND GO TO READ BASE NORMALIZATION"

#Read base normalization of decontaminated reads
echo "########## Normalizing all host nucDNA reads to $(awk -v tar=$target_size "BEGIN {printf \"%.1f\n\", tar / 1000000000}") Gbases"
echo "Paired files: ${sample}.notCombined_1.fastq ${sample}.notCombined_2.fastq"
echo "Single file:  ${sample}_SingleEndQualFiltered.fastq"

echo "calculating sample-base targets..."
p_base=$(cat ${sample}.notCombined_?.fastq | paste - - - - | cut -f2 | tr -d '\n' | wc -c); echo "$p_base bases in paired files"
s_base=$(cat ${sample}_SingleEndQualFiltered.fastq | paste - - - - | cut -f2 | tr -d '\n' | wc -c); echo "$s_base bases in single file"
tot_base=$(awk -v s=$s_base -v p=$p_base "BEGIN {print s + p}"); echo "total $tot_base bases"
paired_target=$(awk -v tot=$tot_base -v p=$p_base -v tar=$target_size "BEGIN {printf \"%.2f\n\", p * tar / tot / 1000000000}")
single_target=$(awk -v tot=$tot_base -v s=$s_base -v tar=$target_size "BEGIN {printf \"%.2f\n\", s * tar / tot / 1000000000}")
echo "target bases for paried reads: $paired_target Gbases"
echo "target bases for single reads: $single_target Gbases"
echo
echo "sampling reads..."   #This utilises 'reformat.sh' from bbmap tools.
echo "paired"
reformat.sh in=${sample}.notCombined_1.fastq in2=${sample}.notCombined_2.fastq\
 out=${sample}.notCombined_norm_1.fastq out2=${sample}.notCombined_norm_2.fastq samplebasestarget=${paired_target}g sampleseed=1
echo "single"
reformat.sh in=${sample}_SingleEndQualFiltered.fastq\
 out=${sample}_SingleEndQualFiltered_norm.fastq samplebasestarget=${single_target}g sampleseed=1
rm -f ${sample}.notCombined_1.fastq ${sample}.notCombined_2.fastq ${sample}_SingleEndQualFiltered.fastq
echo "Normalization (capping) complete"
echo

echo "########## all processes complete"
echo "final read files:"
ls -l #####*NoBacCont*fastq
echo "# compressing and moving results"
mkdir -p $outdir/logs
pigz ${sample}.notCombined_norm_1.fastq ${sample}.notCombined_norm_2.fastq ${sample}_SingleEndQualFiltered_norm.fastq
mv *fastq.gz $outdir 		#final clean reads are moved.
mv bbsplit_*fq.gz $outdir 	#bbsplit reads for symbionts and mitochondria are also moved. Symbiont reads were kept for another study (optional).
rm -r /export/home/a-e/amateos/phyluce_test_2/mapping_tests/scratch/tmp.$SLURM_JOB_ID/*    #delete all remaining files in tmp folder.

echo 
echo "# job ID: $SLURM_JOB_ID finished."
date




