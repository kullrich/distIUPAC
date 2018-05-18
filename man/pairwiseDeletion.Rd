\name{pairwiseDeletion}
\alias{pairwiseDeletion}
\title{pairwiseDeletion}
\usage{
pairwiseDeletion(myvector)
}
\arguments{
\item{myvector}{\code{as.character(DNAStringSet)}}
}
\description{
This function returns the number of used sites for pairwise distance calculations obtained via \link[ape]{dist.dna}.
}
\examples{
data("woodmouse", , package = "ape")
pairwiseDeletion(as.character(dnabin2dnastring(woodmouse)))
}
\author{
Kristian K Ullrich
}
\seealso{
\link[distIUPAC]{distIUPAC}, \link[distIUPAC]{distIUPACmatrix}, \link[ape]{dist.dna}
}
