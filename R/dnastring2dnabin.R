#' @title dnastring2dnabin
#' @name dnastring2dnabin
#' @description This function transforms a \code{DNAStringSet} object from the
#' \code{Biostrings} package into an \code{DNAbin} class of the
#' \code{ape} package.
#' @import Biostrings
#' @importFrom ape as.character.DNAbin as.DNAbin
#' @param dnastring \code{DNAStringSet}
#' @examples
#' data("MySequences", package="distIUPAC")
#' alignment.bin<-dnastring2dnabin(MySequences)
#' @seealso \link[Biostrings]{DNAStringSet}, \link[ape]{as.DNAbin}
#' @export dnastring2dnabin
#' @author Kristian K Ullrich
dnastring2dnabin<-function(dnastring){
    alignment.nb<-length(dnastring)
    alignment.nam<-names(dnastring)
    alignment.seq<-tolower(as.character(dnastring))
    alignment.com<-NA
    alignment<-list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
    names(alignment)<-c("nb", "nam", "seq", "com")
    attr(alignment, "class")<-"alignment"
    alignment.bin<-as.DNAbin(alignment)
    return(alignment.bin)
}
