#' @title dnastring2apealg
#' @name dnastring2apealg
#' @description This function transforms a \code{DNAStringSet} object from the
#' \code{Biostrings} package into an \code{alignment} class of the
#' \code{ape} package.
#' @import Biostrings
#' @importFrom ape as.character.DNAbin as.DNAbin
#' @param dnastring \code{DNAStringSet}
#' @examples
#' data("MySequences", package="distIUPAC")
#' alignment<-dnastring2apealg(MySequences)
#' alignment.bin<-ape::as.DNAbin(alignment)
#' @seealso \link[Biostrings]{DNAStringSet}, \link[ape]{as.alignment}
#' @export dnastring2apealg
#' @author Kristian K Ullrich
dnastring2apealg<-function(dnastring){
    alignment.nb<-length(dnastring)
    alignment.nam<-names(dnastring)
    alignment.seq<-tolower(as.character(dnastring))
    alignment.com<-NA
    alignment<-list(alignment.nb, alignment.nam, alignment.seq, alignment.com)
    names(alignment)<-c("nb", "nam", "seq", "com")
    attr(alignment, "class")<-"alignment"
    return(alignment)
}
