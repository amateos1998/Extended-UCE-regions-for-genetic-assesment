#!/bin/bash

export R_LIBS="/export/home/a-e/amateos/phyluce_test_2/tools/Rlibs"

#This script converts .dist files from ngsDist into MDS matrix files, using the ngsTools script "getMDS.R". 

for infile in $(ls /your/path/results/12_ngs_dist_region/*.dist); do
 echo "Infile: $infile"
outfile=$(echo "$infile" | sed 's|12_ngs_dist_region/|12_ngs_dist_region/mds_dist/mds_dist|' | sed 's/\.dist$/.mds/')
 n_ind=$(head $infile -n2 | tail -n1)
 tail -n+3 $infile | head -n $n_ind | Rscript --vanilla --slave /your/path/ngsTools/Scripts/getMDS.R \
     --no_header --data_symm -n 4 -m "MDS" -o $outfile
 echo
done
