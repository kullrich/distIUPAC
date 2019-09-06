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
#' @param x.pos population X positions [default: NULL]
#' @param min.ind minimum number of individuals without gaps ("-", "+", ".")
#' or without missing sites ("N"), set to size of population ("length(x.pos)")
#' to mask global deletion sites [default: 0]
#' @param wlen sliding window length [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param threads number of parallel threads [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#'   [default: TRUE]
#' @seealso \code{\link[distIUPAC]{triSites}},
#' \code{\link[distIUPAC]{biSites}}
#' @export iupacString
#' @author Kristian K Ullrich
iupacString<-function(dna, name="iupacString", wlen=25000, start.by=1,
  end.by=NULL, threads=1, pB=TRUE){
    options(scipen=22)
    if(length(dna)!=2){stop("dna needs to be of length 2")}
    if(is.null(end.by)){end.by<-unique(width(dna))}
    if(start.by>unique(width(dna))){stop("start.by needs to be equal or
      smaller than dna length")}
    if(end.by>unique(width(dna))){stop("end.by needs to be equal or
      smallerthan dna length")}
    tmp.sw<-swgen(wlen=wlen, wjump=wlen, start.by=start.by, end.by=end.by)
    if(pB){
        pb<-txtProgressBar(min=0, max=ncol(tmp.sw), initial=0, style=3)
    }
    j<-NULL
    registerDoMC(threads)
    OUT<-foreach(j=seq(from=1, to=ncol(tmp.sw)), .combine=c) %dopar% {
        START<-tmp.sw[1, j][[1]]
        END<-tmp.sw[2, j][[1]]
        tmp.seq<-subseq(dna, START, END)
        iString<-rcpp_iupacString_ab(as.character(tmp.seq[1]),
          as.character(tmp.seq[2]), END-START+1, name)
        if(pB){
            setTxtProgressBar(pb, j)
        }
        list(iString)
    }
    if(pB){
        setTxtProgressBar(pb, ncol(tmp.sw))
        close(pb)    
    }
    OUT<-DNAStringSet(paste(OUT, collapse=""))
    names(OUT)<-name
    return(OUT)
}
