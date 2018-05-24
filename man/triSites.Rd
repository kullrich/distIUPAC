% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/triSites.R
\name{triSites}
\alias{triSites}
\title{triSites}
\usage{
triSites(dna, x.pos = NULL, wlen = 25000, threads = 1, pB = TRUE)
}
\arguments{
\item{dna}{\code{DNAStringSet}}

\item{x.pos}{population X positions}

\item{wlen}{sliding window length}

\item{threads}{number of parallel threads}

\item{pB}{specifies if progress should be shown as a progress bar}
}
\description{
This function returns tri-allelic site positions given a dna object
\code{DNAStringSet}, also with IUPAC code.
}
\examples{
data("MySequences", package = "distIUPAC")
#consider all sequences
MySequences.triSites<-triSites(MySequences)
as.matrix(MySequences)[,head(MySequences.triSites)]
#consider only a subset of all sequences
CAS.pos<-5:34
CAS.triSites<-triSites(MySequences, x.pos = CAS.pos)
as.matrix(MySequences[CAS.pos])[,head(CAS.triSites)]
}
\author{
Kristian K Ullrich
}