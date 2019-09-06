\name{rcpp_distIUPAC}
\alias{rcpp_distIUPAC}
\title{rcpp_distIUPAC}
\usage{
rcpp_distIUPAC(dnavector, ncores=1)
}
\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}

\item{ncores}{number of parallel cores to process pairwise distance
calculation [default: 1]}
}
\description{
This function calculates pairwise \code{IUPAC} distances on \code{Biostrings}
\code{DNAStringSet} and \code{BStringSet} objects in a parallel fashion.
}
\examples{
data("MySequences", package="distIUPAC")
rcpp_distIUPAC(as.character(MySequences[1:10]), ncores=2)
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{distIUPAC}},
\code{\link[distIUPAC]{rcpp_distIUPAC_ab}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix}},
\code{\link[ape]{dist.dna}}
}
