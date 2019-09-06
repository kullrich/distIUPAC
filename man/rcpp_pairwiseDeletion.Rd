\name{rcpp_pairwiseDeletion}
\alias{rcpp_pairwiseDeletion}
\title{rcpp_pairwiseDeletion}
\usage{
rcpp_pairwiseDeletion(dnavector, ncores=1)
}
\arguments{
\item{dnavector}{\code{as.character(DNAStringSet)}}

\item{ncores}{number of parallel cores to process pairwise deletion
calculation [default: 1]}
}
\description{
This function returns the number of used sites for pairwise distance
 calculations obtained via \code{\link[ape]{dist.dna}}.
}
\examples{
data("woodmouse", package="ape")
rcpp_pairwiseDeletion(as.character(dnabin2dnastring(woodmouse)), ncores=2)
}
\author{
Kristian K Ullrich
}
\seealso{
\code{\link[distIUPAC]{pairwiseDeletion}},
\code{\link[distIUPAC]{distIUPAC}},
\code{\link[distIUPAC]{distIUPACmatrix}},
\code{\link[ape]{dist.dna}}
}
