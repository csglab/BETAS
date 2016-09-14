#!/bin/bash

echo "Source metadata file: $1"
echo "Source count file: $2"
echo "Sample to be extracted: $3"
echo "New metadata file: $4"
echo "New count file: $5"

cat $1 | cut -f1-3 | cat - $2 | awk -v sampleLabel="$3" -f ./src/_awk/ExtractSample.awk > $5
head -n 1 $1 | cut -f1-3 > $4
cat $1 | awk -v sampleLabel="$3" '$2==sampleLabel' >> $4
echo -e "All.exon\tAll\texonic" >> $4
echo -e "All.intron\tAll\tintronic" >> $4
