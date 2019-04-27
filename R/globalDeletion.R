#' @title globalDeletion
#' @name globalDeletion
#' @description This function returns a \code{DNAStringSet} reduced by all sites containing any gaps ("-", "+", ".") or missing ("N") sites
#' @import Biostrings
#' @import ape
#' @param dna \code{DNAStringSet}
#' @examples
#' data("MySequences", package = "distIUPAC")
#' mySubSequence <- globalDeletion(MySequences)
#' @seealso \link[distIUPAC]{pairwiseDeletion}, \link[distIUPAC]{biSites}, \link[distIUPAC]{triSites}
#' @export globalDeletion
#' @author Kristian K Ullrich
globalDeletion<-function(dna){
  cM<-consensusMatrix(dna)
  globalDeletionSites<-which(apply(cM,2,function(x) sum(x[15:18])>=1))
  if(length(globalDeletionSites)==0){
    return(dna)
  }
  return(dnabin2dnastring(as.DNAbin(DNAMultipleAlignment(dna))[,-globalDeletionSites]))
}