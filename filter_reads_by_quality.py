# script to filter out entire reads if a specified proportion of bases fall below a specified quality

# usage: python3 filter_reads_by_quality.py -i [input_file] -o [output_file] -q [quality_cutoff] -p [percent]
# will output only reads that have at least the given percent of bases at or above the given quality

import argparse
import sys
import gzip

# set up arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help = 'input file name', dest = 'infile')
parser.add_argument('-o', '--output', help = 'output file name', dest = 'outfile')
parser.add_argument('-b', '--base', help = 'ascii base offset (64 or 33)', dest = 'base', type=int, choices=[33,64], default=64)
parser.add_argument('-q', '--quality', help = 'quality score cutoff', dest = 'qual_cut', type=int)
parser.add_argument('-p', '--percent', help = 'percent of read that must make cutoff to keep read', dest = 'percent', type=float)

args = parser.parse_args()

# can use gzipped or uncompressed fastq files
if args.infile.endswith('.gz'):
        infile = gzip.open(args.infile, mode = 'rt')
else:
        infile = open(args.infile)
outfile = open(args.outfile, 'w')

# function to get quality score based on ascii
def get_qual(char):
    qual = ord(char)-args.base
    return(qual)

# function to get sequences of barcodes and reads
def section(infile):
    # read in lines four at a time
    read = [infile.readline(), infile.readline(), infile.readline(), infile.readline()]
    read = [i.strip() for i in read] # remove whitespace
    # return
    if read[3] == "":
        return None
    return read

#print(section(infile))

# get quality score for each base (on last of 4 lines)
while True:
    read = section(infile)
    if read == None:
        break
#    print('\n',  list(read[3]))
    quals = [get_qual(x) for x in list(read[3])]

    # only write reads to output if the proportion of based above the specified quality
    # is higher than the specified proportion
    prop = sum(n >= args.qual_cut for n in quals)/len(quals) # proportion bases meeting cutoff
    if prop*100 >= args.percent:
        output = "\n".join(read) + "\n"
#        print(output)
        outfile.write(output)

infile.close()
outfile.close()
