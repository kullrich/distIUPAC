\name{distIUPAC}
\alias{distIUPAC}
\title{distIUPAC}
\usage{
distIUPAC(myvector)
}
\arguments{
\item{myvector}{\code{as.character(DNAStringSet)}}
}
\description{
This function calculates IUPAC distances on \code{Biostrings} \code{DNAStringSet} and \code{BStringSet} objects.
}
\examples{
data("MySequences", package = "distIUPAC")
distIUPAC(as.character(MySequences[1:10]))
}
\author{
Kristian K Ullrich
}
\seealso{
\link[distIUPAC]{distIUPACmatrix}, \link[ape]{dist.dna}
}
