#!/bin/bash
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_slurm_outplaceholder_id_array[$SLURM_ARRAY_TASK_ID].txt
#SBATCH --error=placeholder_slurm_outplaceholder_id_array[$SLURM_ARRAY_TASK_ID].txt
#SBATCH --job-name rrr_make_loom_files
#SBATCH --array=0-placeholder_max_array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --partition=himem,general
#SBATCH --time=12:00:00
#SBATCH --wait

ml SAMtools/1.15

ml load miniforge3/24.3.0
eval "$(conda shell.bash hook)"

#get arrays
bam_array=(tmp_bams/*bam)
cell_array=(tmp_bcs/*.tsv)

#get current bam file and barcode file from arrays
bam_file=${bam_array[$SLURM_ARRAY_TASK_ID]}
cell_file=${cell_array[$SLURM_ARRAY_TASK_ID]}

echo bam file is $bam_file and cell file is $cell_file

#get sampleid from bam file
#MAKE THIS SAMPLEID placeholder
sampleid=${bam_file:9:5}

#Get arguments from sbatch template
loom_dir=placeholder_loom_dir
env_path=placeholder_env_path
gtf_file=placeholder_gtf_file

#activate environment
conda activate $env_path

#Run velocyto
velocyto run \
    ${bam_file} \
    -b ${cell_file} \
    -o ${loom_dir} \
    -e ${sampleid} \
    ${gtf_file}