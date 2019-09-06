\name{rcpp_distIUPACmatrix_ab}
\alias{rcpp_distIUPACmatrix_ab}
\title{rcpp_distIUPACmatrix_ab}
\usage{
rcpp_distIUPACmatrix_ab(a, b, scoreMatrix, nsites)
}
\arguments{
\item{a}{\code{as.character(DNAString)}}

\item{b}{\code{as.character(DNAString)}}

\item{scoreMatrix}{\code{scoreMatrix}}

\item{nsites}{number of sites}
}
\description{
This function calculates pairwise \code{IUPAC} distance between
two \code{character} objects (\code{a}, \code{b})of length \code{nsites} using
a specific \code{scoreMatrix} and returns the number of used sites.
}
\examples{
##load sequence data
data("MySequences", package="distIUPAC")
##load scoring matrix
MyScoreMatrix<-scoreMatrix()
rcpp_distIUPACmatrix_ab(as.character(MySequences[1]), as.character(MySequences[2]),
scoreMatrix=MyScoreMatrix, nsites=nchar(as.character(MySequences[1]))[[1]])
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{scoreMatrix}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix}},
\code{\link[ape]{dist.dna}}
}
