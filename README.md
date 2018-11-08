# RNAseq_scripts

Set of scripts for processing RNAseq data files.

## convert_qseq_to_fastq.py:
Converts qseq sequence files to fastq files.

Usage: python3 convert_qseq_to_fastq -i [input qseq file] -o [output fastq file] [options]

Written in python 3.

Assumes qseq files are formatted like those from UCLA's BSCRC sequencing core: tab-separated file, columns are:
instrument name, run ID, lane number, tile number, x coordinate, y coordinate, index, end number, read, quality scores, and filter (1=passes, 0 = fails).

Can take gzipped or uncompressed files as input.

Can filter out reads that fail the Illumina quality filter (default), or keep them in the output (option --nofilter).

## demultiplexer.py:
Demultiplexes fastq files - reads with the same barcode sequence are output to the same fastq files.

Usage: python3 demultiplexer.py ['set_A' or 'set_B'] ['paired' or 'single'] [path to directory with files]

This works with paired-end or single-end reads, and gzipped or uncompressed files.

It will identify Illumia TruSeq LT set A or set B indices, depending on which you specify. (Info about barcode sequences here: http://web.archive.org/web/20160327094802/http://support.illumina.com:80/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseqsampleprep/truseq-library-prep-pooling-guide-15042173-01.pdf)


The script assumes names of fastq files are formatted like this:

s_&ast;_1_&ast;.fastq&ast; file with read 1 sequence

s_&ast;_2_&ast;.fastq&ast; file with barcodes

s_&ast;_3_&ast;.fastq&ast; file with read 2 sequence (if paired-end)


The script will include barcodes with 1 mismatched base in the demultiplexed file for that barcode. Because barcode sequences from UCLA's BSCRC sequencing core include 7 bases, the script will match the full 7-base barcode, the barcode with a "." in the first position instead of the correct character (this seems to be fairly common), or the 7-base barcode with one of the bases incorrect. For example, all of these barcode sequences will match index 25:

ACTGATA (the actual barcode)

.CTGATA (the barcode with a "." in the first position)

ATTGATA (the barcode with one base mismatch in the second position)


The script will output files named after the barcode sequence: ACTGATA_read1.fastq (and ACTGATA_read2.fastq if paired-end), as well as a file for the unmatched reads: NOMATCH_read1.fastq and NOMATCH_read2.fastq.
