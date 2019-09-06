\name{distIUPACmatrix}
\alias{distIUPACmatrix}
\title{distIUPACmatrix}
\usage{
distIUPACmatrix(dnavector, scoreMatrix)
}\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}

\item{scoreMatrix}{\code{scoreMatrix}}
}
\description{
This function calculates scaled pairwise \code{IUPAC} distances on
\code{Biostrings} \code{DNAStringSet} and \code{BStringSet} objects using a
specific \code{scoreMatrix} in a sequential fashion.
}
\examples{
##load sequence data
data("MySequences", package="distIUPAC")
##load scoring matrix
MyScoreMatrix<-scoreMatrix()
distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
##change scoring matrix
MyScoreMatrix["A", "R"]<-10.0
##calculate IUPAC distance with altered scoring
distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{scoreMatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix_ab}}
\code{\link[ape]{dist.dna}}
}
