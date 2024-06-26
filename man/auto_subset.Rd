% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scSeurat.R
\name{auto_subset}
\alias{auto_subset}
\title{Use median and standard deviation to subset a Seurat object based on specific features}
\usage{
auto_subset(
  sobject,
  sd_down = 1,
  sd_up = 2,
  make_plots = TRUE,
  features = c("nCount_RNA", "nFeature_RNA"),
  sample_name = NULL
)
}
\arguments{
\item{sobject}{Seurat object}

\item{sd_down}{Number of standard deviations below the median to subset}

\item{sd_up}{Number of standard deviations above the median to subset}

\item{make_plots}{Whether to make plots of the features before and after subsetting}

\item{features}{Vector of features to use for subsetting}

\item{sample_name}{Name of sample to use in plot titles}
}
\value{
A Seurat object
}
\description{
Use median and standard deviation to subset a Seurat object based on specific features
}
\examples{
\dontrun{
filtered <- auto_subset(SeuratObject::pbmc_small)
}
}
