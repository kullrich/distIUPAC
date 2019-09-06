\name{distIUPAC}
\alias{distIUPAC}
\title{distIUPAC}
\usage{
distIUPAC(dnavector)
}
\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}
}
\description{
This function calculates pairwise \code{IUPAC} distances on \code{Biostrings}
\code{DNAStringSet} and \code{BStringSet} objects
in a sequential fashion.
}
\examples{
data("MySequences", package="distIUPAC")
distIUPAC(as.character(MySequences[1:10]))
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{rcpp_distIUPAC}},
\code{\link[distIUPAC]{rcpp_distIUPAC_ab}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix}},
\code{\link[ape]{dist.dna}}
}
