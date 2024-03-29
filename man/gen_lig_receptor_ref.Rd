% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scIntercellular.R
\name{gen_lig_receptor_ref}
\alias{gen_lig_receptor_ref}
\title{Get reference data for find_ligand() function}
\usage{
gen_lig_receptor_ref(
  lig_tar_matrix = "https://zenodo.org/record/3260758/files/ligand_target_matrix.rds",
  lig_rec_network = "https://zenodo.org/record/3260758/files/lr_network.rds",
  weighted_network = "https://zenodo.org/record/3260758/files/weighted_networks.rds"
)
}
\arguments{
\item{lig_tar_matrix}{Ligand target matrix R data file location}

\item{lig_rec_network}{Ligand receptor network R data file location}

\item{weighted_network}{Ligand receptor weighted network R data file location}
}
\value{
None
}
\description{
Get reference data for find_ligand() function
}
\keyword{internal}
