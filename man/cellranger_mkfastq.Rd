% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProcessData.R
\name{cellranger_mkfastq}
\alias{cellranger_mkfastq}
\title{Run cellranger mkfastq}
\usage{
cellranger_mkfastq(
  sample_info,
  email = "",
  tar_folders,
  bcl_folder = "/home/gdrobertslab/lab/BCLs",
  fastq_folder = "/home/gdrobertslab/lab/FASTQs",
  slurm_out = paste(getwd(), "/slurmOut_mkfastq-\%j.out", sep = "")
)
}
\arguments{
\item{sample_info}{File containing sample info (see details)}

\item{email}{Email for Slurm notifications}

\item{tar_folders}{Character vector containing folders extracted from tar}

\item{bcl_folder}{Path to write BCL files}

\item{fastq_folder}{Path to write fastq files}

\item{slurm_out}{Location to write out slurm out files}
}
\description{
Run cellranger mkfastq
}
\examples{
\dontrun{
need example
}
}
