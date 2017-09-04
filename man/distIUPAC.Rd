\name{distIUPAC}
\alias{distIUPAC}
\title{distIUPAC}
\usage{
distIUPAC()
}
\description{
}
\examples{
data("MySequences", package = "distDNA")
distIUPAC(as.character(MySequences$dna))
distIUPAC(as.character(MySequences$iupac))
}
\author{
Kristian K Ullrich
}
\seealso{
\link[distIUPAC]{distIUPACmatrix}, \link[ape]{dist.dna}
}
