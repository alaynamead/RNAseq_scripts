
# script to convert qseq to fastq files
# usage: -i [infile] -o [outfile] [options]

import argparse
import sys
import gzip

# set up arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'input (qseq) file name', dest = 'infile')
parser.add_argument('-o', '--output', help = 'output (fastq) file name', dest = 'outfile')
filter_parser = parser.add_mutually_exclusive_group(required = False)
filter_parser.add_argument('--filter', help='will not include reads that did not pass Illumina quality filter in output (default)', dest = 'filter', action = 'store_true')
filter_parser.add_argument('--no_filter', help='will include all reads in output', dest = 'filter', action = 'store_false')
parser.set_defaults(filter = True)

args = parser.parse_args()

# can use gzipped or normal qseq files
if args.infile.endswith('.gz'):
    infile = gzip.open(args.infile, mode = 'rt')
else:
    infile = open(args.infile)
outfile = open(args.outfile, 'w')

# make function to convert a line of a seq file to fastq format
# line must be converted to list first
def convert(read):
    # define name of read
    readname = ':'.join(read[0:8] + list(read[10]))
    output = '@{}\n{}\n+\n{}\n'.format(readname, read[8], read[9])
    return output

# convert each line (read) from qseq file to 4-line fastq format
for line in infile:

    # check if lines should be filtered or not
    if args.filter:
        # convert each column to element in a list
        line_list = line.split()
        # only write lines if Illumina filter was passed (last column = 1)
        if(int(line_list[10]) == 1):
            outfile.write(convert(line_list))

    # if not filtering, don't need to check if filter was passed
    else:
        line_list = line.split()
        outfile.write(convert(line_list))

infile.close()
outfile.close()
