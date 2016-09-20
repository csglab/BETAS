#!/bin/bash

metadata="./metadata.txt"
countfile="./ENCSR648QFY_PCBP2_shRNA.counts.txt"
betas="../../bin/BETAS"

outdir="./betas.out"
mkdir -p $outdir
outfile=$outdir"/ENCSR648QFY_PCBP2_shRNA.BETAS"

$betas -meta $metadata -counts $countfile -out $outfile > $outfile.cout.log 2>$outfile.cerr.log
