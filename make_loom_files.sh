#!/bin/bash
ml SAMtools/1.15
ml load Miniconda3/4.9.2

eval "$(conda shell.bash hook)"

#Get location of rrrSingleCellUtils and create environment there

#conda prefix

# Check that all user input ID's have BAM files and check that the user
# wants to continue if not

# Decide what gtf file to use based on species

OPTSTRING=":x:o:s:"

while getopts ${OPTSTRING} opt; do
    case ${opt} in
        x)
            echo "Option -x was triggered, argument: ${OPTARG}"
            samples=${OPTARG}
            ;;
        o)
            echo "Option -o was triggered, argument: ${OPTARG}"
            output_path=${OPTARG}
            ;;
        s)
            echo "Option -s was triggered, argument: ${OPTARG}"
            species=${OPTARG}
            ;;
        b)
            echo "Option -b was triggered, argument: ${OPTARG}"
            barcode_path=${OPTARG}
            ;;
        :)
            echo "Option -${OPTARG} requires an argument"
            exit 1
            ;;
        ?)
            echo "INVALID OPTION: -${OPTARG}"
            exit 1
            ;;
    esac
done

#This checks if any of the arguments are missing and exits if so
if [[ -z "${samples}" || -z "${output_path}" || -z "${species}" ]]; then
    echo -e "Missing mandatory arguments\n-x is a string with sample ID's separated by a space\n-o is the directory in which you want your loom files saved to\n-s is the species from which your sample came from. Must be either \"human\" or \"mouse\"\nExiting..."
    exit 1
fi

if [[ -z "${barcode_path}" ]]; then
    echo "No barcode path provided. Using all barcodes from 10x output"

#Split samples text string to array
###TESTING IT OUT
IFS=' ' read -a samples <<< "$samples"

#Set variables to make running velocyto easier
#path to barcodes (will change later)
bc_path_0="/home/gdrobertslab/lab/Analysis/MattGust/projects/Roberts_Lab/f420_single_nuclei/misc/"
bc_path_1="_bcs.tsv"

#Set path to bam files
bam_path_0="/home/gdrobertslab/lab/Counts_2/"
bam_path_1="/possorted_genome_bam.bam"

#Create appropriate path to gtf file based on input species
echo "species is ${species}"
if [[ ${species} == "human" ]]; then
    gtf_path="/home/gdrobertslab/lab/GenRef/10x-hg38/genes/genes.gtf.gz"
elif [[ ${species} == "mouse" ]]; then
    gtf_path="/home/gdrobertslab/lab/GenRef/10x-mm10/genes/genes.gtf.gz"
else
    echo "Unknown species. Species parameter must be either \"human\" or \"mouse\". Exiting..."
    exit 1
fi

#loop through array of samples
for id in ${samples[@]} 
do
    old_bam="${bam_path_0}${id}${bam_path_1}"
    #make temporary folder for bams and then copy bam file into it
    mkdir tmp_bam
    cp $old_bam tmp_bam/$id
    args=`echo "bam_files/$id -b ${bc_path_0}${id}${bc_path_1} -o $output_path $gtf_path"`
    velocyto run $args
    rm -r tmp_bam
done