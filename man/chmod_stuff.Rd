% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getProcessData.R
\name{chmod_stuff}
\alias{chmod_stuff}
\title{Set folder permissions to read-only.}
\usage{
chmod_stuff(folder, code = "444", record_log = TRUE)
}
\arguments{
\item{folder}{The path to the folder.}

\item{code}{The code to use for chmod.}

\item{record_log}{Record failures when processing data}
}
\value{
None
}
\description{
This function takes a folder path as input and sets the permissions of the
 folder to read-only.
}
\keyword{internal}
