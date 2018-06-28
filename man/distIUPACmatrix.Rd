\name{distIUPACmatrix}
\alias{distIUPACmatrix}
\title{distIUPACmatrix}
\usage{
distIUPACmatrix(myvector, scoreMatrix)
}\arguments{
\item{myvector}{\code{as.character(DNAStringSet)}}

\item{scoreMatrix}{\code{scoreMatrix}}
}
\description{
This function calculates IUPAC distances on \code{Biostrings} \code{DNAStringSet} and \code{BStringSet} objects and a \code{scoreMatrix}.
}
\examples{
data("MySequences", package = "distIUPAC")
MyScoreMatrix<-scoreMatrix()
distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
MyScoreMatrix["A","R"] <- 10.0
distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
MyScoreMatrix["A","R"] <- 0.5
distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
}
\author{
Kristian K Ullrich
}
\seealso{
\link[distIUPAC]{scoreMatrix}, \link[ape]{dist.dna}
}
