import subprocess
from multiprocessing import Pool
import argparse
import re
import sys
import pandas as pd
import numpy as np
import itertools
np.seterr(invalid='ignore')

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--vcf',
                    type=str,
                    default='test.bcf',
                    help='VCF file with multiple samples as columns')
parser.add_argument('--out_base',
                    '-o',
                    type = str,
                    default = 'out_dist',
                    help='output tsv file name')
parser.add_argument('--verbose',
                    action='store_true',
                    help='print out extra information')
parser.add_argument('--trim_path',
                    action='store_true',
                    help='trim path from sample names')
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
dist_key_dict = {'0/0': 0,
                 '0/1': 1,
                 '1/0': 1,
                 '1/1': 2,
                 './.': np.nan}

################################################################################
### Code

########
### functions

# Calculate the distance between two samples using dist_key_dict and add to dict

# Process a batch of vcf lines to create a list of two pandas dataframes (dist and count)
def process_line(vcf_line):
    vcf_list = [re.sub(':.+', '', x) for x in vcf_line.strip().split('\t')]
    # Kick out any SNP with multiple alt calls
    if ',' not in vcf_list[4]:
        # replace all instances of dist_key_dict keys with the values
        vcf_list = [dist_key_dict.get(ele, np.nan) for ele in vcf_list[9:]]
        dist = abs(np.array(vcf_list) - np.array(vcf_list).reshape(-1, 1))
        # Need to multiply by 2 since we're comparing two alleles at each locus
        count = (np.isnan(dist) == False).astype('float64') * 2

        return(np.array([np.nan_to_num(dist, 0), count]))
    else:
        return(empty_array)

def process_lines(chr):
    cmd_list = ['bcftools', 'view', '-Hr', chr, args.vcf]
    output = sum(map(process_line,
                     subprocess.check_output(cmd_list).\
                        decode("utf-8").\
                        splitlines()))
    if args.verbose:
        print(chr, file = sys.stderr)

    if type(output) != int:
        return(output)
    else:
        return(empty_array)

### Check in vcf/bcf index exists


########
# Read in the VCF header
cmd_list = ['bcftools', 'view', '-h', args.vcf]
vcf_data = subprocess.\
            check_output(cmd_list).\
            decode("utf-8").\
            splitlines()

batch_size = 1000000
# Read in the file until the header
for line in vcf_data:
    if line.startswith(r'##contig=<ID='):
        line = re.sub(r'##contig=<ID=', '', line.strip())
        chr_len = int(re.sub(r'.+length=', '', line).strip(">"))
        starts = list(range(1, chr_len, batch_size))
        stops = [x - 1 for x in starts]
        stops.append(chr_len)
        stops = stops[1:]
        chr_name = re.sub(r',.+', '', line)
        chr_list.extend([i + ":" + str(j) + "-" + str(k) for i, j, k in \
            zip(itertools.repeat(chr_name), starts, stops)])
    if line.startswith('#CHROM'):
        header = line.strip().split('\t')
        header = header[9:]
        if args.trim_path:
            header = [re.sub('.+/', '', x) for x in header]
        break # exit the loop since the next line will start the data

empty_array =\
    np.array([np.zeros(len(header)**2).reshape(len(header),
                                               len(header)),
              np.zeros(len(header)**2).reshape(len(header),
                                               len(header))])

# For the test of the file:
# Make list of lists of lines to pass to pool.map

with Pool(processes = args.processes) as pool:
    output = sum(pool.imap(process_lines, chr_list))


# Calculate the average distance
ave_dist = output[0, :] / output[1, :]
# Make it a square matrix
ave_dist.shape = (len(header), len(header))
# Make it a dataframe with the header as the index and column names
ave_dist = pd.DataFrame(ave_dist,
                        index = header,
                        columns = header)
# write ave_dist_df to file
ave_dist.to_csv(args.out_base + '.tsv', sep = '\t', na_rep = 'NA')

tot_dist_df = output[0, :]
tot_dist_df.shape = (len(header), len(header))
tot_dist_df = pd.DataFrame(tot_dist_df,
                           index = header,
                           columns = header)
tot_dist_df.to_csv(args.out_base + '_tot_dist.tsv', sep = '\t', na_rep = 'NA')

tot_count_df = output[1, :]
tot_count_df.shape = (len(header), len(header))
tot_count_df = pd.DataFrame(tot_count_df,
                            index = header,
                            columns = header)
tot_count_df.to_csv(args.out_base + '_tot_count.tsv', sep = '\t', na_rep = 'NA')
