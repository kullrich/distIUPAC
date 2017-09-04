\name{distIUPACmatrix}
\alias{distIUPACmatrix}
\title{distIUPACmatrix}
\usage{
distIUPACmatrix()
}
\description{
}
\examples{
data("MySequences", package = "distIUPAC")
MyScoreMatrix<-scoreMatrix()
distIUPACmatrix(as.character(MySequences$iupac), MyScoreMatrix)
MyScoreMatrix["A","R"] <- 10.0
distIUPACmatrix(as.character(MySequences$iupac), MyScoreMatrix)
MyScoreMatrix["A","R"] <- 0.5
distIUPACmatrix(as.character(MySequences$iupac), MyScoreMatrix)
}
\author{
Kristian K Ullrich
}
\seealso{
\link[distIUPAC]{scoreMatrix}, \link[ape]{dist.dna}
}
