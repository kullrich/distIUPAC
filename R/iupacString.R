#' @title iupacString
#' @name iupacString
#' @description This function returns a \code{IUPAC} sequence based on a
#' \code{Biostrings} objects
#' @import Biostrings
#' @import foreach
#' @import doMC
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet} of length 2 [mandatory]
#' @param x.pos population X positions [default: 1]
#' @param y.pos population X positions [default: 2]
#' @param wlen sliding window length [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param threads number of parallel threads [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#'   [default: TRUE]
#' @seealso \code{\link[distIUPAC]{triSites}},
#' \code{\link[distIUPAC]{biSites}}
#' @examples
#' data("MySequences", package = "distIUPAC")
#' iupacString(MySequences, x.pos=5, y.pos=6)
#' @export iupacString
#' @author Kristian K Ullrich
iupacString<-function(dna, x.pos=1, y.pos=2, name="iupacString",
  wlen=25000, start.by=1, end.by=NULL, threads=1, pB=FALSE){
    options(scipen=22)
    dna_<-dna[c(x.pos, y.pos)]
    if(length(dna_)!=2){stop("dna needs to be of length 2")}
    OUT<-tmpSEQsw(dna_,
      FUN=function(x) {
          x.len<-unique(width(x))
          list(
            rcpp_iupacString_ab(as.character(x[1]), as.character(x[2]),
              x.len, "iupacString")
          )
      }, chr.name="chr", wlen=wlen, wjump=wlen,
      start.by=start.by, end.by=end.by, wtype="bp",
      global.deletion=FALSE, threads=threads, pB=pB)
    OUT.seq<-DNAStringSet(paste0(OUT[,4],collapse = ""))
    names(OUT.seq)<-name
    return(OUT.seq)
}
