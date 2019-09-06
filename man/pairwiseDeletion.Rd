\name{pairwiseDeletion}
\alias{pairwiseDeletion}
\title{pairwiseDeletion}
\usage{
pairwiseDeletion(dnavector)
}
\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}
}
\description{
This function returns the number of used sites for pairwise distance
 calculations obtained via \code{\link[ape]{dist.dna}}.
}
\examples{
data("woodmouse", package="ape")
pairwiseDeletion(as.character(dnabin2dnastring(woodmouse)))
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{rcpp_pairwiseDeletion}},
\code{\link[distIUPAC]{distIUPAC}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[ape]{dist.dna}}
}
