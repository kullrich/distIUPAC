% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnastring2dnabin.R
\name{dnastring2dnabin}
\alias{dnastring2dnabin}
\title{dnastring2dnabin}
\usage{
dnastring2dnabin(dnastring)
}
\arguments{
\item{dnastring}{\code{DNAStringSet}}
}
\description{
This function transforms a \code{DNAStringSet} object from the
\code{Biostrings} package into an \code{DNAbin} class of the
\code{ape} package.
}
\examples{
data("MySequences", package="distIUPAC")
alignment.bin<-dnastring2dnabin(MySequences)
}
\seealso{
\link[Biostrings]{DNAStringSet}, \link[ape]{as.DNAbin}
}
\author{
Kristian K Ullrich
}
