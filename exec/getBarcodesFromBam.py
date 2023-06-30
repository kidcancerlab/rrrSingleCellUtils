import argparse
import sys
import re
# import pandas as pd
# import numpy as np
import fcntl
import subprocess
import itertools
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Get sam entries for each cell barcode.')
parser.add_argument('--cells',
                    '-c',
                    type=str,
                    default='F420_cells_10cells.txt',
                    help='list of cell barcodes with no header')
parser.add_argument('--bam',
                    '-b',
                    type=str,
                    default='F420_100k.bam',
                    help='bam file output from cellranger')
parser.add_argument('--out_base',
                    '-o',
                    type = str,
                    default = 'cellbam_',
                    help='output bam file name')
parser.add_argument('--sam_batch_n',
                    '-n',
                    type = int,
                    default=10000000,
                    help='number of sam entries to process at a time')
parser.add_argument('--verbose',
                    action='store_true',
                    help='print out extra information')
parser.add_argument('--processes',
                    '-p',
                    type = int,
                    default=1,
                    help='number of processes to use for parallel processing')

args = parser.parse_args()

################################################################################
### Global variables
chr_list = []
header = []
read_names = []
all_reads = {}

################################################################################
### Code

########
### functions

#########
# Process a batch of bam lines to create a dictionary of sam entries for each barcode
def process_lines(chrom):
    lines_dict = {}
    umi_dict = {}
    cmd_list = ['samtools', 'view', args.bam, chrom]
    # check if any reads in range provided by chrom
    if len(subprocess.check_output(cmd_list).\
                    decode("utf-8").\
                    splitlines()) != 0:
        for line in subprocess.check_output(cmd_list).\
                        decode("utf-8").\
                        splitlines():
            if re.compile(r'CB:Z:([ATGC]+-1)\t').search(line) and \
                   re.compile(r'CB:Z:([ATGC]+-1)\t').search(line).group(1) in cell_barcodes:
                cell_barcode = re.compile(r'CB:Z:([ATGC]+-1)\t').search(line).group(1)
                umi_chr_loc = re.compile(r'UB:Z:([ATGC]+)\t').\
                                 search(line).\
                                 group(1) + \
                    '_' + \
                    line.split('\t')[2] + \
                    '_' + \
                    line.split('\t')[3]

                # Add barcode to list inside lines_dict
                if cell_barcode not in lines_dict:
                    lines_dict[cell_barcode] = []

                # Don't add line if umi_chr_loc + cell_barcode already in umi_dict
                # This gets rid of PCR duplicates
                if umi_chr_loc + cell_barcode not in umi_dict:
                    lines_dict[cell_barcode].append(line)
                    umi_dict[umi_chr_loc + cell_barcode] = 1

    def print_dict(barcode):
        nonlocal lines_dict
        with open(args.out_base + barcode + ".sam", "a") as g:
            fcntl.flock(g, fcntl.LOCK_EX)
            g.writelines('\n'.join(lines_dict.get(barcode)) + '\n')
            fcntl.flock(g, fcntl.LOCK_UN)
        return

    # Print out the dictionary for each barcode
    [print_dict(barcode) for barcode in list(lines_dict.keys())]
    # print_dict(barcode)
    if args.verbose:
        print(chrom + " is done.", file = sys.stderr)

    return

def write_header(barcode):
    open(args.out_base + barcode + ".sam", "w").writelines('\n'.join(bam_header) + '\n')

########
# Read in cell barcodes
# This assumes there is no header in barcode file
barcode_file = open(args.cells, 'r')
cell_barcodes = [x.strip() for x in barcode_file.readlines()]
# check if cell_barcodes is empty
if not cell_barcodes:
    print("Error: cell_barcodes is empty", file = sys.stderr)
    sys.exit(1)

########
# Read in the bam header
cmd_list = ['samtools', 'view', '-H', args.bam]
bam_header = subprocess.\
            check_output(cmd_list).\
            decode("utf-8").\
            splitlines()

# batch_size = 1000000
# Process header to
for line in bam_header:
    if line.startswith(r'@SQ'):
        chr_len = int(re.sub(r'.+\tLN:', '', line.strip()))
        starts = list(range(1, chr_len, args.sam_batch_n))
        stops = [x - 1 for x in starts]
        stops.append(chr_len)
        stops = stops[1:]
        chr_name = re.compile(r'SN:(.+?)\t').search(line).group(1)
        chr_list.extend([i + ":" + str(j) + "-" + str(k) for i, j, k in \
            zip(itertools.repeat(chr_name), starts, stops)])

# Print out length of chr_list
if args.verbose:
    print(str(len(chr_list)) + " chunks to process", file = sys.stderr)

# Print out the header into each output sam file
with Pool(processes = args.processes) as pool:
    output = pool.map(write_header, cell_barcodes)

# For the test of the file:
# Make list of lists of lines to pass to pool.map

with Pool(processes = args.processes) as pool:
    output = pool.map(process_lines, chr_list)

output = list(output)
print(len(output))
