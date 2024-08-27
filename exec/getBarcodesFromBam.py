import pysam
import argparse
import sys
import os

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
parser.add_argument('--verbose',
                    action='store_true',
                    help='print out extra information')

args = parser.parse_args()

################################################################################
### Global variables
label_dict = {}
umi_dict = {}
################################################################################
### Code

########
### functions

def add_to_label_dict(x):
    cell, label = x.split('\t')
    label_dict[cell] = label

def open_bam_outs_from_labels(labels, bam_template):
    bam_outs = {}
    for label in labels:
        # check if the output bam file already exists
        bam_file = args.out_base + '_' + label + '.bam'
        if os.path.exists(bam_file):
            print("Error: bam file already exists", file = sys.stderr)
            sys.exit(1)
        else:
            bam_outs[label] = \
                pysam.AlignmentFile(bam_file, 'wb', template = bam_template)
    return bam_outs


# For a single bam line, get the CB:Z and UB:Z tags, then if the CB:Z tag is in the
# label_dict as a key, write the line to the appropriate bam file if the UB:Z
# tag has not been seen before (not in the umi_dict dictionary)
def process_line(line, bam_outs):
    # check if CB:Z tag is present
    if line.has_tag('CB') and line.has_tag('UB'):
        # get CB:Z and UB:Z tags
        cb = line.get_tag('CB')
        ub = line.get_tag('UB')
        molecule = cb + ub + str(line.reference_name) + str(line.reference_start)
        # check if CB:Z tag is in label_dict
        if cb in label_dict:
            # check if UB:Z tag is in umi_dict
            if molecule not in umi_dict:
                # add UB:Z tag to umi_dict
                umi_dict[molecule] = 1
                # write line to bam file
                bam_outs[label_dict[cb]].write(line)
        return(1)
    else:
        return(0)

# make main function
def main():
    in_bam = pysam.AlignmentFile(args.bam, 'rb')
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

    bam_outs = open_bam_outs_from_labels(all_labels, in_bam)

    # Loop through the bam file and use process_line on each line
    # # Use multiprocessing to parallelize
    # pool = Pool(processes = args.processes)
    # pool.map(process_line, in_bam, itertools.repeat(bam_outs))
    # pool.close()
    counter = 0
    good_lines = 0
    for line in in_bam:
        good_lines += process_line(line, bam_outs)
        if counter % 1000000 == 0 and args.verbose:
            print('{:,d}'.format(counter) + \
                    ' reads processed. ' + \
                    '{:,d}'.format(good_lines) + \
                    " kept.",
                  file=sys.stderr)
        counter += 1

if __name__ == '__main__':
    main()