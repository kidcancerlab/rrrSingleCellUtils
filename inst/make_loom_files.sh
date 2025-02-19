#!/bin/bash
#SBATCH --account=placeholder_account
#SBATCH --output=placeholder_slurm_out%j.txt
#SBATCH --error=placeholder_slurm_error%j.txt
#SBATCH --job-name rrr_make_loom_files
#SBATCH --array=0-placeholder_max_array%50
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=240G
#SBATCH --partition=himem,general
#SBATCH --time=16:00:00
#SBATCH --wait

ml SAMtools/1.15

ml load miniforge3/24.3.0
eval "$(conda shell.bash hook)"

#get output directory
out_dir=placeholder_out_dir
    
#get arrays
id_array=(placeholder_id_array)

#get sampleid from id array
sampleid=${id_array[$SLURM_ARRAY_TASK_ID]}

#get current barcode file from out_dir
cell_file=${out_dir}${sampleid}/tmp_bcs.tsv


#make bam file from out_dir and sampleid
bam_file=${out_dir}${sampleid}/tmp_bam.bam

echo bam file is $bam_file and cell file is $cell_file

#rename output file with sample id
echo This is job ID ${SLURM_JOBID}
echo Replacing ${SLURM_JOBID}.txt with ${sampleid}.txt

mv placeholder_slurm_out${SLURM_JOBID}.txt \
    ${out_dir}${sampleid}/slurm_out.txt

mv placeholder_slurm_error${SLURM_JOBID}.txt \
    ${out_dir}${sampleid}/slurm_error.txt

#Get arguments from sbatch template
loom_dir=${out_dir}${sampleid}
env_path=placeholder_env_path
gtf_file=${out_dir}${sampleid}/genes.gtf.gz

gunzip ${gtf_file} ${out_dir}${sampleid}/genes.gtf

gtf_file=${out_dir}${sampleid}/genes.gtf

#activate environment
conda activate $env_path

#Run velocyto
velocyto run \
    ${bam_file} \
    -b ${cell_file} \
    -o ${loom_dir} \
    -e ${sampleid} \
    ${gtf_file}

echo "Just ran velocyto"

# rm -rf ${out_dir}${sampleid}_bams/
rm -rf ${out_dir}${sampleid}/genes*
rm -rf ${out_dir}${sampleid}/*.bam
echo "Just deleted bams. Job's done"

sleep 60


mv placeholder_slurm_out%j.txt ${out_dir}/${sampleid}/job_output_%j.txt