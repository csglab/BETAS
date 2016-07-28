# BETAS
##Bayesian Estimation of Transcript Abundance and Stability

BETAS implements Bayesian statistics to estimate the distribution of abundance of exonic and intronic fragments of each gene from the observed exonic and intronic read counts. It uses Markov Chain Monte Carlo simulations with Metropolis-Hastings algorithm to estimate the parameters of distribution.

### Installation

From the command line, go to the root folder of BETAS, and run `make`. This will compile the codes and will create the BETAS executable in `./bin/`.

### Running BETAS


#### Input files

For running BETAS, you need the following files:

* **Metadata file**. The metadata file is a tab-delimited table with each row corresponding to either the exonic or intronic reads of one sample. Lines that start with `#` will be treated as comments, and will not be read. There should be at least three columns in the metadata file, in the same order as follows:
	1. **Read set name**. This should be a unique label. For example, you can label the intronic reads from sample #1 as `Sample1.intronic`, and the exonic reads from sample #1 as `Sample1.exonic`.
	2. **Sample name**. For example, both the above read sets that come from sample #1 should be labeled as `Sample1`.
	3. **Read type**. Should be either `intronic` or `exonic`.
* **Read count table**. In the read count table, each row corresponds to one gene, and each column corresponds to one of the read sets specified in the metadata file. The column headers should exactly match the read set names in the metadata file. See below regarding how to create the read count table from htseq-count output files.

##### Creating read count table from htseq-count files

In order to create the read count table, all htseq-count files from intronic and exonic read sets of all samples need to be merged. I suggest including an additional column (column 4) in the metadata file, corresponding to the path of the htseq-count file that contains the counts for each read set. An example table is shown below:

| # Read Set     | Label   | Read Type | HTSeq-Count File                        |
| -------------- | ------- | --------- | --------------------------------------- |
| Sample1.exon   | Sample1 | exonic    | ./htseq/Sample1.htseqCount.exonic.tab   |
| Sample1.intron | Sample1 | intronic  | ./htseq/Sample1.htseqCount.intronic.tab |
| Sample2.exon   | Sample2 | exonic    | ./htseq/Sample2.htseqCount.exonic.tab   |
| Sample2.intron | Sample2 | intronic  | ./htseq/Sample2.htseqCount.intronic.tab |

Then, use a shell script like this to create the read count table:

```bash

metadata="<metadat_file.txt>"
countFile="<readcount_file.txt>"

# the perl script join.pl can merge multiple two-column tab-delimited files based on IDs in the first column
multijoin="<BETAS_folder>/src/_perl/join.pl"

for i in `cat $metadata | sed 1d | cut -f1`
do
	# get the htseq-count file name that corresponds to this sample
	htseqFile=`cat $metadata | grep -w $i | cut -f4`
	
	# remove the last 5 lines, which contain the counts of the reads that do not map to a gene
	# NOTE: this does not work on all platforms, since some versions of "head" may not accept negative numbers
	head -n -5 $htseqFile > $i
done

# merge the files, and include only the genes that appear in all files
$multijoin -c `cat $metadata | sed 1d | cut -f1` > $countFile

# remove the temporary count files
rm -f `cat $metadata | sed 1d | cut -f1`
```


#### Usage

To run BETAS, use the following command:

```bash
./bin/BETAS -meta \<metadata_file.txt\> -counts \<readcount_file.txt\> -out \<output_file.out.txt\>
```

BETAS has other parameters for tuning its behaviour, but they are experimental at this moment, and should not be modified.

#### Output

BETAS creates an output table, in which each row corresponds to one gene, and each column corresponds to one of the following variables for each sample:

* `fi_average`: The estimated log of abundance of intronic fragments.
* `fi_stdev`: The standard deviation of the above estimate (confidence intervals).
* `fe_average`: The estimated log of abundance of exonic fragments.
* `fe_stdev`: The standard deviation of the above estimate (confidence intervals).
* `dfi_average`: The estimated log ratio of abundance of intronic fragments relative to other samples (Δintron).
* `dfi_stdev`: The standard deviation of the above estimate (confidence intervals).
* `dfe_average`: The estimated log ratio of abundance of exonic fragments relative to other samples (Δexon).
* `dfe_stdev`: The standard deviation of the above estimate (confidence intervals).
* `ddf_average`The estimated Δexon–Δintron (change in stability).
* `ddf_stdev`: The standard deviation of the above estimate (confidence intervals).

The above columns are repeated in the output file for each sample, with the sample name indicated on the first line of the file on top of each column.

#### Example

An example dataset is provided at `./examples/mouse.GM_2D`. Simply go to this folder, and run `bash run.BETAS.sh`. It should automatically create the read count table from htseq-count files that are in `./examples/mouse.GM_2D/htseq`, and run BETAS. The output will be written to `./examples/mouse.GM_2D/betas.out`. Two log files will also be created; the file `mouse.GM_2D.BETAS.cerr.log` will contain the acceptance rate for the MCMC simulation (useful for debugging), and the file `mouse.GM_2D.BETAS.cout.log` will contain the main messages that are printed by BETAS. The main output file of BETAS with the estimated abundances will be `mouse.GM_2D.BETAS.out.txt`.
