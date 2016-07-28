#!/bin/bash

metadata="./metadata.txt"
multijoin="../../src/_perl/join.pl"
betas="../../bin/BETAS"

### First, convert the HTSeq-Count output to the read count table

countfile="./mouse.GM_2D.counts.txt"

for i in `cat $metadata | sed 1d | cut -f1`
do
    # get the htseq-count file name that corresponds to this sample
    htseqFile=`cat $metadata | grep -w $i | cut -f4`

    # remove the last 5 lines, which contain the counts of the reads that do not map to a gene
    # NOTE: this does not work on all platforms, since some versions of "head" may not accept negative numbers
    head -n -5 $htseqFile > $i
done

# merge the files, and include only the genes that appear in all files
$multijoin -c `cat $metadata | sed 1d | cut -f1` > $countfile

# remove the temporary count files
rm -f `cat $metadata | sed 1d | cut -f1`

### Now, run BETAS

outdir="./betas.out"
mkdir -p $outdir
outfile=$outdir"/mouse.GM_2D.BETAS"

$betas -meta $metadata -counts $countfile -out $outfile > $outfile.cout.log 2>$outfile.cerr.log

rm -f $countfile
