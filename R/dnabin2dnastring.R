#' @title dnabin2dnastring
#' @name dnabin2dnastring
#' @description This function transforms a \code{DNAbin} object from the \code{ape} package into an \code{DNAStringSet} class of the \code{Biostrings} package.
#' @import Biostrings
#' @import ape
#' @param dna \code{DNAbin}
#' @examples
#' data("woodmouse", package = "ape")
#' alignment <- dnabin2dnastring(woodmouse)
#' @seealso \link[Biostrings]{DNAStringSet}, \link[ape]{DNAbin}
#' @export dnabin2dnastring
#' @author Kristian K Ullrich
dnabin2dnastring<-function(dna){
  alignment<-DNAStringSet(apply(as.character(dna),1,function(x) paste0(x,collapse="")))
  return(alignment)
}
