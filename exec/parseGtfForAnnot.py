import argparse
import gzip


parser = argparse.ArgumentParser(description='Process gtf file to get just the bits I need.')
parser.add_argument('--gtf',
                    type=str,
                    help='gtf file')

args = parser.parse_args()

################################################################################
### Global variables
wanted_fields = ['gene_biotype', 'transcript_id', 'gene_id', 'gene_name']
obs_fields = set()
################################################################################
### Code

########
# Read in the gtf file and print out the bits I need
if args.gtf.endswith('.gz'):
    gtf_file = gzip.open(args.gtf, mode = 'r')
else:
    gtf_file = open(args.gtf, 'r')

# print header
print('gene_biotype', 'tx_name', 'gene_id', 'gene_name', sep = '\t')

for line in gtf_file:
    line = line.decode('utf-8').strip()
    output = []
    if not line.startswith('#'):
        col_list = line.split('\t')
        # putting these fields into a dict since I don't know if all the ones I want will be there
        gtf_dict = dict()
        info_list = col_list[8].strip(';').split('; ')
        for field in info_list:
            key, value = field.split(' ', maxsplit = 1)
            gtf_dict[key] = value.strip('"')
        for field in wanted_fields:
            if field in gtf_dict:
                output.append(gtf_dict[field])
            else:
                output.append('NA')
        # only keep novel combinations
        if ('\t'.join(output) not in obs_fields) and output[1] != 'NA':
            print('\t'.join(output))
            obs_fields.add('\t'.join(output))

gtf_file.close()