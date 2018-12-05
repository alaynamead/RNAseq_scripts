# script to demultiplex fastq files
# written by Adam Fontenot and Alayna Mead

# usage: 
# python3 demultiplexer.py ['set_A' or 'set_B'] ['paired' or 'single'] [path to directory with files]

import sys
import glob
import gzip

# arguments
bar_set = sys.argv[1] # which barcode set to use ('set_A' or 'set_B')
paired = sys.argv[2] # paired or single end reads ('paired' or 'single')
path = sys.argv[3] # path to directory with read files

if bar_set not in ('set_A', 'set_B'):
    sys.exit('first argument (adapter set) must be set_A or set_B')

if paired not in ('paired', 'single'):
    sys.exit('second argument (read type) must be paired or single')

# get files into a list
# are files paired or single end reads?
if paired == 'paired':
    seq_locs_1 = glob.glob(path + '/s_*_1_*.fastq*') # file with read 1 sequence
    bar_locs = glob.glob(path + '/s_*_2_*.fastq*') # file with barcodes
    seq_locs_2 = glob.glob(path + '/s_*_3_*.fastq*') # file with read 2 sequence

    # sort list of files alphabetically, so order is the same
    seq_locs_1.sort()
    bar_locs.sort()
    seq_locs_2.sort()

    if not len(seq_locs_1) == len(seq_locs_2) == len(bar_locs):
        sys.exit('There are not the same number of read 1, read 2, and barcode files')

elif paired == 'single':
    seq_locs_1 = glob.glob(path + '/s_*_1_*.fastq*') # file with read 1 sequence
    bar_locs = glob.glob(path + '/s_*_2_*.fastq*') # file with barcodes
    # sort list of files
    seq_locs_1.sort()
    bar_locs.sort()

    if not len(seq_locs_1) == len(bar_locs):
        sys.exit('There are not the same number of sequence and barcode files')


# which barcode set to use?
# Illumina set A
if bar_set == 'set_A':
    barcodes = ["CGATGTA", "TGACCAA", "ACAGTGA", "GCCAATA", "CAGATCA", "CTTGTAA",
               "AGTCAAC", "AGTTCCG", "ATGTCAG", "CCGTCCC", "GTCCGCA", "GTGAAAC"]
    barcode_names = ["index_02", "index_04", "index_05", "index_06", "index_07", "index_12",
                    "index_13", "index_14", "index_15", "index_16", "index_18", "index_19"]

# Illumina set B
if bar_set == 'set_B':
    barcodes = ["ATCACGA", "TTAGGCA", "ACTTGAA", "GATCAGA", "TAGCTTA", "GGCTACA",
               "GTGGCCT", "GTTTCGG", "CGTACGT", "GAGTGGA", "ACTGATA", "ATTCCTT"]
    barcode_names = ["index_01", "index_03", "index_08", "index_09", "index_10", "index_11",
                    "index_20", "index_21", "index_22", "index_23", "index_25", "index_27"]


letters = set("ACGT")

# this sets up a dictionary with all acceptable "matches" to a barcode
# by switching out each base with another
# allows 1 base mismatch to barcode, or barcode starting with "." instead of base
barcode_lookup = {}
for barcode in barcodes:
    barcode_lookup[barcode] = barcode
    barcode_lookup['.'+barcode[1:]] = barcode
    for i in range(len(barcode)):
        for l in letters - set(barcode[i]):
            barcode_lookup[barcode[:i] + l + barcode[i+1:]] = barcode

# function to get sequences of barcodes and reads
def section(seq_file_1, bar_file, seq_file_2=""):
    
    seqs_1 = [seq_file_1.readline(),seq_file_1.readline(),seq_file_1.readline(),seq_file_1.readline()] # read lines
    seqs_1 = [i.strip() for i in seqs_1] # remove whitespace
    bars = [bar_file.readline(),bar_file.readline(),bar_file.readline(),bar_file.readline()]
    bars = [i.strip() for i in bars]
    
    if seq_file_2 != "": # do same if there is a read 2 sequence file (for paired end reads)
        seqs_2 = [seq_file_2.readline(),seq_file_2.readline(),seq_file_2.readline(),seq_file_2.readline()]
        seqs_2 = [i.strip() for i in seqs_2]
    
    # return
    if seq_file_2 == "": # if there is no read 2 file (for single end reads)
        if seqs_1[3] == "":
              return None, None
        return seqs_1, bars
        
    if seq_file_2 != "": # paired end
        if seqs_1[3] == "":
            return None, None, None
        return seqs_1, bars, seqs_2

# function to return actual barcode from any barcode sequence matching it (ie barcode with a mismatch)
def match(barcode):
    try:
        return barcode_lookup[barcode]
    except:
        return None

# write each read to file based on which barcode it matches

# single end:
if paired == 'single':

    output_files = {}
    for seq_loc, bar_loc in zip(seq_locs_1, bar_locs):
        # can read gz or uncompressed files
        # assumes sequence and barcode files for each tile are both gz or not
        if seq_loc[-3:] == '.gz':
            seq_file = gzip.open(seq_loc, mode = 'rt')
            bar_file = gzip.open(bar_loc, mode = 'rt')
        else:
            seq_file = open(seq_loc, 'r')
            bar_file = open(bar_loc, 'r')
        
        no_match_file = open("NOMATCH.fastq", "a") # file for reads that don't match a barcode

        # make new file for each barcode
        for x in range(0,len(barcodes)):
            filename = barcode_names[x] + ".fastq" 
            output_files[barcodes[x]] = open(filename, "a")
        
        while True:
            seqs, bars = section(seq_file, bar_file)
            if seqs == None:
                break
            barcode = bars[1] # actual barcode is on 2nd line of 4 lines
            barcode = match(barcode)
            if barcode == None:
                 output = "\n".join(seqs) + "\n"
                 no_match_file.write(output)
            else:
                #seqs[0] = seqs[0] + " forward " + barcode
                output = "\n".join(seqs) + "\n"
                output_files[barcode].write(output)        
                
        # close files
        seq_file.close()
        bar_file.close()
        no_match_file.close()
        for output_file in output_files:
            output_files[output_file].close()
        
        print("finished", seq_loc)

elif paired == 'paired':
    
    output_files_read1 = {}
    output_files_read2 = {}
   
    for seq_loc_1, seq_loc_2, bar_loc, in zip(seq_locs_1, seq_locs_2, bar_locs):
        # can read gz or uncompressed files
        # assumes sequence and barcode files for each tile are all gz or not
        if seq_loc_1[-3:] == '.gz':
            seq_file_1 = gzip.open(seq_loc_1, mode = 'rt')
            seq_file_2 = gzip.open(seq_loc_2, mode = 'rt')
            bar_file = gzip.open(bar_loc, mode = 'rt')
        else:
            seq_file_1 = open(seq_loc_1, 'r')
            seq_file_2 = open(seq_loc_2, 'r')
            bar_file = open(bar_loc, 'r')
   
        no_match_file_1 = open("NOMATCH_read1.fastq", "a") # file for reads that don't match a barcode
        no_match_file_2 = open("NOMATCH_read2.fastq", "a")
        
        # make 2 new files - read 1 and read 2 for each barcode 
        for x in range(0, len(barcodes)):
            filename1 = barcode_names[x] + "_read1.fastq"
            filename2 = barcode_names[x] + "_read2.fastq"
            output_files_read1[barcodes[x]] = open(filename1, "a")
            output_files_read2[barcodes[x]] = open(filename2, "a")

        while True:
            seqs_1, bars, seqs_2 = section(seq_file_1, bar_file, seq_file_2)
            if seqs_1 == None:
                break
            barcode = bars[1] # actual barcode is on 2nd line of 4 lines
            barcode = match(barcode)
            if barcode == None:
                 output_1 = "\n".join(seqs_1) + "\n"
                 no_match_file_1.write(output_1)
                 output_2 = "\n".join(seqs_2) + "\n"
                 no_match_file_2.write(output_2)

            else:
                output_1 = "\n".join(seqs_1) + "\n"
                output_files_read1[barcode].write(output_1)
                output_2 = "\n".join(seqs_2) + "\n"
                output_files_read2[barcode].write(output_2)

        # close files
        seq_file_1.close()
        seq_file_2.close()
        bar_file.close()
        no_match_file_1.close()
        no_match_file_2.close()
        for output_file in output_files_read1:
            output_files_read1[output_file].close()
        for output_file in output_files_read2:
            output_files_read2[output_file].close()
        
        print("finished", seq_loc_1)

