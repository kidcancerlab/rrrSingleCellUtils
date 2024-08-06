#!/bin/bash
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_slurm_out%j.txt
#SBATCH --error=placeholder_slurm_out%j.txt
#SBATCH --job-name rrr_make_loom_files
#SBATCH --array=0-placeholder_max_array%10
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
id_array=(placeholder_id_array)

#get current bam file and barcode file from arrays
bam_file=${bam_array[$SLURM_ARRAY_TASK_ID]}
cell_file=${cell_array[$SLURM_ARRAY_TASK_ID]}

echo bam file is $bam_file and cell file is $cell_file

#get sampleid from id array
sampleid=${id_array[$SLURM_ARRAY_TASK_ID]}

#rename output file with sample id
echo This is job ID ${SLURM_JOBID}
echo Replacing ${SLURM_JOBID}.txt with ${sampleid}.txt

mv placeholder_slurm_out${SLURM_JOBID}.txt \
    placeholder_slurm_out${sampleid}.txt

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

rm ${bam_file}