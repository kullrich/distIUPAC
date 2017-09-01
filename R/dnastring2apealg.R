#' @title dnastring2aperalg
#' @name dnastring2apealg
#' @description This function transforms a \code{DNAStringSet} object from the \code{Biostrings} package into an \code{alignment} class of the \code{ape} package.
#' @import Biostrings
#' @param dna \code{DNAStringSet}
#' @examples
#' data("MySequences", package = "distDNA")
#' alignment <- dnastring2apealg(MySequences$dna)
#' alignment.bin <- as.DNAbin(alignment)
#' @seealso \link[Biostrings]{DNAStringSet}, \link[ape]{as.alignment}
#' @export dnastring2apealg
#' @author Kristian K Ullrich
dnastring2apealg<-function(dna){
  alignment.nb<-length(dna)
  alignment.nam<-names(dna)
  alignment.seq<-tolower(as.character(dna))
  alignment.com<-NA
  alignment<-list(alignment.nb,alignment.nam,alignment.seq,alignment.com)
  names(alignment)<-c("nb","nam","seq","com")
  attr(alignment,"class")<-"alignment"
  return(alignment)
}
