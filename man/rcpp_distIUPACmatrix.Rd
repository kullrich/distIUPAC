\name{rcpp_distIUPACmatrix}
\alias{rcpp_distIUPACmatrix}
\title{rcpp_distIUPACmatrix}
\usage{
rcpp_distIUPACmatrix(dnavector, scoreMatrix, ncores=1)
}\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}

\item{scoreMatrix}{\code{scoreMatrix}}

\item{ncores}{number of parallel cores to process pairwise distance
calculation [default: 1]}
}
\description{
This function calculates scaled pairwise \code{IUPAC} distances on
\code{Biostrings} \code{DNAStringSet} and \code{BStringSet} objects using a 
specific \code{scoreMatrix} in a parallel fashion.
}
\examples{
##load sequence data
data("MySequences", package="distIUPAC")
##load scoring matrix
MyScoreMatrix<-scoreMatrix()
rcpp_distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix, ncores=2)
##change scoring matrix
MyScoreMatrix["A","R"] <- 10.0
##calculate IUPAC distance with altered scoring
rcpp_distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix, ncores=2)
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{scoreMatrix}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[distIUPAC]{rcpp_distIUPACmatrix_ab}},
\code{\link[ape]{dist.dna}}
}
