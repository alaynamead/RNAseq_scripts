# RNAseq_scripts

Set of scripts for processing RNAseq data files.

## convert_qseq_to_fastq.py:
Usage: python3 convert_qseq_to_fastq -i [input qseq file] -o [output fastq file] [options]

Written in python 3.

Assumes qseq files are formatted like those from UCLA's BSCRC sequencing core: tab-separated file, columns are:
instrument name, run ID, lane number, tile number, x coordinate, y coordinate, index, end number, read, quality scores, and filter (1=passes, 0 = fails).

Can take gzipped or uncompressed files as input.

Can filter out reads that fail the Illumina quality filter (default), or keep them in the output (option --nofilter).
