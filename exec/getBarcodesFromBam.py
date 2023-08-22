import argparse
import sys
import re
import fcntl
import subprocess
import time
import itertools
from multiprocessing import Pool
from multiprocessing import Lock

parser = argparse.ArgumentParser(description='Get sam entries for each cell barcode.')
parser.add_argument('--cells',
                    '-c',
                    type=str,
                    default='/home/gdrobertslab/mvc002/analyses/roberts/dev/testSnps/output/two_cancers/two_cancers_clusters_cells_S0149.txt',
                    help='list of cell barcodes with no header')
parser.add_argument('--bam',
                    '-b',
                    type=str,
                    default='/home/gdrobertslab/lab/Counts/S0149/possorted_genome_bam.bam',
                    help='bam file output from cellranger')
parser.add_argument('--out_base',
                    '-o',
                    type = str,
                    default = 'test',
                    help='output bam file name')
parser.add_argument('--sam_batch_n',
                    '-n',
                    type = int,
                    default=100000000,
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
label_dict = {}
chr_list = []
header = []
read_names = []
all_reads = {}
lock = Lock()

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
            line = line.strip()
            if re.compile(r'CB:Z:([ATGC]+-1)').search(line) and \
               re.compile(r'CB:Z:([ATGC]+-1)').search(line).group(1) in label_dict.keys() and \
               re.compile(r'UB:Z:([ATGC]+)').search(line):
                cell_barcode = re.compile(r'CB:Z:([ATGC]+-1)\t').search(line).group(1)
                umi = re.compile(r'UB:Z:([ATGC]+)').\
                                 search(line).\
                                 group(1)
                chr = line.split('\t')[2]
                pos = line.split('\t')[3]
                molecule = \
                    umi + \
                    chr + \
                    pos + \
                    cell_barcode

                cell_label = label_dict.get(cell_barcode)

                # Add barcode to list inside lines_dict
                if cell_label not in lines_dict:
                    lines_dict[cell_label] = []

                # Don't add line if molecule already in umi_dict
                # This gets rid of PCR duplicates
                if molecule not in umi_dict:
                    lines_dict[cell_label].append(line)
                    umi_dict[molecule] = 1

    def print_dict(label):
        nonlocal lines_dict
        while True:
            try:
                with open(args.out_base + label + ".sam", "a") as filehandle:
                    fcntl.flock(filehandle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    filehandle.writelines('\n'.join(lines_dict.get(label)) + '\n')
                    fcntl.flock(filehandle, fcntl.LOCK_UN)
                    filehandle.close()
                    break
            except IOError:
                print("Error: " + label + " is locked", file = sys.stderr)
                time.sleep(1)
        return

    # Print out the dictionary for each barcode
    with lock:
        junk = [print_dict(label) for label in list(lines_dict.keys())]

    if args.verbose:
        print(chrom + " is done.", file = sys.stderr)

    return

def write_header(label):
    sam_file = open(args.out_base + label + ".sam", "w")
    sam_file.writelines('\n'.join(bam_header) + '\n')
    sam_file.close()

def add_to_label_dict(x):
    cell, label = x.split('\t')
    label_dict[cell] = label

########
# Read in cell barcodes
# This assumes there is no header in barcode file
barcode_file = open(args.cells, 'r')
cell_barcodes = [add_to_label_dict(x.strip()) for x in barcode_file.readlines()]
all_labels = list(set(label_dict.values()))
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
    output = pool.map(write_header, all_labels)

# For the test of the file:
# Make list of lists of lines to pass to pool.map

with Pool(processes = args.processes) as pool:
    output = pool.map(process_lines, chr_list)

output = list(output)
print(len(output))
