#' @title globalDeletion
#' @name globalDeletion
#' @description This function returns a \code{DNAStringSet} reduced by all
#' sites containing any gaps ("-", "+", ".") or missing ("N") sites.
#' @import Biostrings
#' @importFrom ape as.character.DNAbin as.DNAbin
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions [default: NULL]
#' @param wlen sliding window length [default: 50000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param threads number of parallel threads [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @examples
#' data("MySequences", package="distIUPAC")
#' mySubSequence<-globalDeletion(MySequences)
#' @seealso \link[distIUPAC]{pairwiseDeletion}, \link[distIUPAC]{biSites},
#' \link[distIUPAC]{triSites}
#' @export globalDeletion
#' @author Kristian K Ullrich
globalDeletion<-function(dna, x.pos=NULL, wlen=50000, start.by=1, end.by=NULL,
  threads=1, pB=FALSE){
    options(scipen=22)
    if(is.null(end.by)){end.by<-unique(width(dna))}
    if(start.by>unique(width(dna))){
        stop("start.by needs to be equal or smaller than dna length")
    }
    if(end.by>unique(width(dna))){
        stop("end.by needs to be equal or smaller than dna length")
    }
    if(is.null(x.pos)){
        x.pos<-seq(1, length(dna))
    }
    dna_<-dna[x.pos]
    x.pos_<-seq(1, length(x.pos))
    tmp.sw<-swgen(wlen=wlen, wjump=wlen, start.by=start.by, end.by=end.by)
    if(pB){
        pb<-txtProgressBar(min=0, max=ncol(tmp.sw), initial=0, style=3)
    }
    j<-NULL
    registerDoMC(threads)
    OUT<-foreach(j=seq(from=1, to=ncol(tmp.sw)), .combine=c) %dopar% {
        START<-tmp.sw[1, j][[1]]
        END<-tmp.sw[2, j][[1]]
        tmp.seq<-subseq(dna_, START, END)
        cM<-consensusMatrix(dna)
        globalDeletionSites<-which(apply(cM, 2, function(x) sum(x[15:18])>=1))
        if(pB){
            setTxtProgressBar(pb, j)
        }
        globalDeletionSites
    }
    if(pB){
        setTxtProgressBar(pb, ncol(tmp.sw))
        close(pb)    
    }
    if(length(OUT)==0){
        return(dna)
    }
    return(dnabin2dnastring(as.DNAbin(DNAMultipleAlignment(dna))[, -OUT]))
}