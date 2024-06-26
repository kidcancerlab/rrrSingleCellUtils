% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProcessData.R
\name{fix_sample_sheet}
\alias{fix_sample_sheet}
\title{Fix sampleSheet.csv for cellranger mkfastq}
\usage{
fix_sample_sheet(
  orig_sample_sheet,
  new_sample_sheet,
  sample_info_sheet,
  record_log = TRUE
)
}
\arguments{
\item{new_sample_sheet}{Location to write new sample sheet}

\item{sample_info_sheet}{Tab delimited sheet with columns containing info to
insert into the sample sheet}

\item{record_log}{Record failures when processing data}

\item{sample_sheet}{SampleSheet.csv from the BCL folder}
}
\description{
Fix sampleSheet.csv for cellranger mkfastq
}
\keyword{internal}
