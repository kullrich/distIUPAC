\name{rcpp_distIUPAC_ab}
\alias{rcpp_distIUPAC_ab}
\title{rcpp_distIUPAC_ab}
\usage{
rcpp_distIUPAC_ab(a, b, nsites)
}
\arguments{
\item{a}{\code{as.character(DNAString)}}

\item{b}{\code{as.character(DNAString)}}

\item{nsites}{number of sites}
}
\description{
This function calculates pairwise \code{IUPAC} distance between
two \code{character} objects (\code{a}, \code{b})of length \code{nsites} and
returns the number of used sites.
}
\examples{
data("MySequences", package="distIUPAC")
rcpp_distIUPAC_ab(as.character(MySequences[1]), as.character(MySequences[2]),
nsites=nchar(as.character(MySequences[1]))[[1]])
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{distIUPAC}},
\code{\link[distIUPAC]{rcpp_distIUPAC}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix_ab}},
\code{\link[ape]{dist.dna}}
}
