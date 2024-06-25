#!/bin/bash
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_slurm_out
#SBATCH --error=placeholder_slurm_out
#SBATCH --job-name rrr_make_loom_files
#SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --partition=himem,general
#SBATCH --time=12:00:00
#SBATCH --wait

echo "abt to load modules"

ml SAMtools/1.15

echo "just loaded samtools"

ml load miniforge3/24.3.0
eval "$(conda shell.bash hook)"

echo "just loaded miniforge"

#Get arguments from sbatch template or whatever
bam_file=placeholder_bam_file
loom_dir=placeholder_loom_dir
cell_file=placeholder_cell_file
env_path=placeholder_env_path
gtf_file=placeholder_gtf_file

echo "just got arguments from sbatch template"

#Check if conda environment exists
if conda info --envs | grep -q ${env_path};
then
    echo "env already exists";
else
    conda create -p ${env_path} -f make_environment.yml
fi

echo "about to activate env"

#activate environment
conda activate $env_path

echo "environment activated, making velocyto call"

#Run velocyto
velocyto run \
    ${bam_file} \
    -b ${cell_file} \
    -o ${loom_dir} \
    ${gtf_file}