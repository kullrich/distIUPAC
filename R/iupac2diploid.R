#' @title iupac2diploid
#' @name iupac2diploid
#' @description This function returns a \code{DNAStringSet} sequence based on a
#' \code{iupacString} diploid \code{Biostrings} object by splitting each diploid
#' sequence into two haploid sequences (no random sampling per site, just the
#' extreme haploid possibilities are used)
#' @import Biostrings
#' @import foreach
#' @import doMC
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet} [mandatory]
#' @seealso \code{\link[distIUPAC]{triSites}},
#' \code{\link[distIUPAC]{biSites}}
#' @examples
#' data("MySequences", package = "distIUPAC")
#' iupac2diploid(MySequences[1])
#' @export iupac2diploid
#' @author Kristian K Ullrich
iupac2diploid<-function(dna){
    options(scipen=22)
    OUT.seq<-DNAStringSet(unlist(lapply(dna, function(x)
      rcpp_iupacString2diploidString(as.character(x), length(x),
        paste0(names(x),"1"), paste0(names(x),"2")))))
    return(OUT.seq)
}
