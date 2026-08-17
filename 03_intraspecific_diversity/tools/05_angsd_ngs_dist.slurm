#!/bin/bash
#SBATCH --job-name=ngsdist
#SBATCH --output=ngsdist_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
#SBATCH --time=9-48:00:00
#SBATCH --partition=cpuq


#angsd_ngs_dist.qsub located in bin
#This script calls SNPs from the indel-realigned sorted BAM files using ANGSD.0.929, GL2 (GATK) and a region filter.
#Resulting genotype likelihoods and posterior probability is used in ngsDist (ngsTools) to calculate distant matrix for phylogenetic tree (optional) and MDS.
#The tree is build by FastME from the distant matrix.
#This script automatically specifies the depth parameters (5-80%tiles) and $minind (50% of samples, rounded up).

reffile=$1    #Oalg_verC3 etc.
sampleinfo=$2 #data/sample_info_180810_310_Oalg.list etc.
bamdir=$3     #results/10b-11_indel_angsd_depth/180810 etc.
outdir=$4     #results/12_ngs_dist etc.
snp_pval=$5   #1e-3 etc.
prior=$6      #1 or 2 (1=observed allele frequency with HWE-assumed, or 2=uniform prior).

####min and max depth can be set automatically with the hashed sections below or manually.
mindepth=250   #####$(sed -e '1,/Global_depth/d' $bamdir/$reffile/*info| awk '$2>0.05 {print $1}'| head -n1) #15 etc.  5 percentile taken from $bamdir/$reffile/*info
maxdepth=450   ######$(sed -e '1,/Global_depth/d' $bamdir/$reffile/*info| awk '$2<0.80 {print $1}'| tail -n1) #75 etc. 80 percentile taken from $bamdir/$reffile/*info

minind=$(cat $sampleinfo | wc -l | awk '{print int($1/2+0.5)}')   #43 etc. 50% of #samples, rounded up integral
maxdepind=10  #Used to filter exome-wide SNP-sites with high sample-depths to exclude repeat-regions in the genome; -1 means no filtering.
hwepval=-1 #Used to monitor SNP-sites extremely departed from HWE; -1 means no filtering on posterior genotype likelihoods.
upperdepth=$(cat $sampleinfo | wc -l | awk '{print $1*10}')   #900 etc. 10X of #samples as the upper limit of global depth-graph
pdir=$(pwd)


####add at the end of angsd command as -rf $regions
regions=/your/path/a_rf.txt    #A list of regions ('chr:') for -rf in ANGSD. Manually made based on 30 samples with ca. 1X coverage to exclude reads from highly repeated regions but accommodate
    # stochasticity. Contigs that have the mean of average coverage <10X are included, but ones that have no mapped reads or all samples had >4X coverage were excluded.
    # Also regions that any of populaitons had no mapped reads were also excluded (common only). 
	# See ../../1_transcriptome_assembly/results/Final_reference_transcriptomes/README.txt for details.


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

#Copy files

mkdir -p /your/path/scratch/tmp.$SLURM_JOB_ID
cd /your/path//scratch/tmp.$SLURM_JOB_ID
echo "####reference files copying#####"
cp /your/path/reference/ahya_20k.fasta $reffile.fasta
echo
echo "Reference: $reffile copied"
echo

echo "Source: $bamdir"
echo "Reference: $reffile"
echo
echo "bam.gz and bai.gz files copying and unzipping..."
for sample in $(cat $sampleinfo | cut -f1); do
 echo $sample
 cp $bamdir/$reffile/$sample*_MinQ20_PairedAndSingle_dedup_sorted_realigned.bam.gz .
done
unpigz *.gz

for sample in $(cat $sampleinfo | cut -f1); do
 cp $bamdir/$reffile/$sample*_MinQ20_PairedAndSingle_dedup_sorted_realigned.bai.gz .
done
unpigz *.gz

ls *.bam > bam_list.txt
n_ind=$(cat bam_list.txt | wc -l)          #Calculate the number of individual samples
cat $sampleinfo | cut -f1 | sort > sample.list
echo "done lists made containing $n_ind samples"
echo



echo "mindepth=$mindepth ; maxdepth=$maxdepth ; minind=$minind will be used in ANGSD"
echo "ANGSD starting..."
date
samtools faidx $reffile.fasta
/your/path/tools/ngsTools/angsd/angsd -P $SLURM_CPUS_PER_TASK -b bam_list.txt -ref $reffile.fasta -out ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}\
            -uniqueOnly 0 -remove_bads 0 -only_proper_pairs 0 -doQsDist 1 -doDepth 1 -maxDepth $upperdepth  -dumpCounts 1\
            -minQ 2 -minInd $minind -setMinDepth $mindepth -setMaxDepth $maxdepth -setMaxDepthInd $maxdepind -doCounts 1\
            -GL 2 -doMajorMinor 1 -doMaf 1 -minMaf 0.01 -SNP_pval $snp_pval\
            -doGeno 8 -doPost $prior -doHWE 1 -minHWEpval $hwepval -rf $regions
Rscript /your/path/tools/ngsTools/Scripts/plotQC.R\
 ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior} 2> /dev/null
echo "ANGSD complete"
NSITES=$(zcat ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}.mafs.gz | tail -n+2 | wc -l)
echo "$NSITES SNP-sites found in ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}.mafs"
date
echo

echo "####activating conda env with gsl2.8 fro ngsDist####"
conda init
conda activate gsl_env
export LD_LIBRARY_PATH=$(conda info --base)/envs/gsl_env/lib:$LD_LIBRARY_PATH

echo "ngsDist starting..."  #phylogenetic tree is optional
echo "single"
/your/path/tools/ngsTools/ngsDist/ngsDist -verbose 1 -geno ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}.geno.gz\
 -probs -n_ind $n_ind -n_sites $NSITES\
 -labels sample.list  -o ${n_ind}Samples_to_${reffile}_MinQ20_SNPpval${snp_pval}_Prior${prior}.dist -n_threads $SLURM_CPUS_PER_TASK
echo "#####Skipping phylogenetic tree######"
echo
echo "##### deactivating conda env with gsl2.8####"
conda deactivate
echo



echo "moving result to:"
echo " $outdir/$reffile"


mkdir -p $outdir/$reffile
mv * $outdir/$reffile
echo
echo
#######rm -r /scratch/tmp.$SLURM_JOB_ID
echo "job ID: $JOB_ID finished."
date
