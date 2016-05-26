#!/bin/bash

metadata="./metadata.txt"
extract="../../src/_awk/ExtractCounts.awk"
multijoin="../../src/_perl/join.pl"
betas="../../bin/BETAS"

### First, convert the HTSeq-Count output to a format suitable for BETAS

countfile="./mouse.GM_2D.counts.txt"

for i in `cat $metadata | sed 1d | cut -f1`
do
	readType=`cat $metadata | grep -w $i | cut -f3 | head -n 1`
	srcFile=`cat $metadata | grep -w $i | cut -f4 | head -n 1`
	
	awk -f $extract -v ReadType=$readType -v FS='\t' $srcFile > $i
done

$multijoin -c `cat $metadata | sed 1d | cut -f1` > $countfile

rm -f `cat $metadata | sed 1d | cut -f1`

### Now, run BETAS

dist=0
learn=1
sampling=2000

outdir="./betas.out"
mkdir -p $outdir
outfile=$outdir"/mouse.GM_2D.BETAS"

$betas -meta $metadata -counts $countfile -dist $dist -learn $learn -samp $sampling -out $outfile > $outfile.cout.log 2>$outfile.cerr.log

rm -f $countfile
