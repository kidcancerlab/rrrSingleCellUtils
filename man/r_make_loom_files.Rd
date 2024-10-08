% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/velocity_analysis_functions.R
\name{r_make_loom_files}
\alias{r_make_loom_files}
\title{Wrapper to Make Loom Files in R}
\usage{
r_make_loom_files(
  sobj,
  sobj_name = NULL,
  out_dir = "loom_output/",
  id_col = NULL,
  species,
  bam_paths,
  cluster_account,
  slurm_base = paste0(getwd(), "/slurmOut"),
  sbatch_base = "sbatch_"
)
}
\arguments{
\item{sobj}{Seurat object you want to run velocity analysis}

\item{sobj_name}{Optional, name of your seurat object for specifying output
when the same sample has a loom created for multiple objects.}

\item{out_dir}{Folder to output loom files and sbatch output to. By default
will create a directory in your working directory called loom_output}

\item{id_col}{Name of metadata column marking what sample a given cell is
from.}

\item{species}{String communicating what species the cells are from.
Value should be either "human" or "mouse".}

\item{bam_paths}{A named list or vector containing the bam files for each
unique sample in id_col. Names should be the sample name and should match the
unique values in sobj[[id_col]]}

\item{cluster_account}{Your Franklin cluster user ID.}

\item{slurm_base}{The directory to write slurm output files to.}

\item{sbatch_base}{The prefix to use with the sbatch job file.}
}
\value{
A loom file for each unique ID present in id_col, output to loom_dir.
}
\description{
Wrapper to Make Loom Files in R
}
\details{
This is a wrapper function around velocyto's "run" command. It will
take your seurat object and generate loom files using only the cells present
in your seurat object. The loom files can then be used to run velocity
analysis using scVelo.
}
