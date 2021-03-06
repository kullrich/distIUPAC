#' @title triSites
#' @name triSites
#' @description This function returns tri-allelic site positions given
#' a dna object \code{DNAStringSet}, also with IUPAC code.
#' @import Biostrings
#' @import foreach
#' @import doMC
#' @importFrom stats as.dist sd setNames
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions [default: NULL]
#' @param min.ind minimum number of individuals without gaps ("-", "+", ".")
#' or without missing sites ("N"), set to size of population ("length(x.pos)")
#' to mask global deletion sites [default: 0]
#' @param wlen sliding window length [default: 25000]
#' @param start.by optional start position[default: 1]
#' @param end.by optional end position [default: NULL]
#' @param threads number of parallel threads [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @seealso \code{\link[distIUPAC]{biSites}},
#' \code{\link[distIUPAC]{globalDeletion}}
#' @examples
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' 
#' ##consider all sequences
#' MySequences.triSites<-triSites(MySequences)
#' as.matrix(MySequences)[,head(MySequences.triSites$triPOS)]
#' 
#' ##consider only a subset of all sequences
#' CAS.pos<-5:34
#' CAS.triSites<-triSites(MySequences, x.pos=CAS.pos)
#' as.matrix(MySequences[CAS.pos])[, head(CAS.triSites$triPOS)]
#' 
#' ##consider only sites were 15 individuals have no gaps or missing sites
#' CAS.triSites.minInd.15<-triSites(MySequences, x.pos=CAS.pos, min.ind=15)
#' as.matrix(MySequences[CAS.pos])[, head(CAS.triSites.minInd.15$triPOS)]
#' as.matrix(MySequences[CAS.pos])[, head(CAS.triSites.minInd.15$triPOSmasked)]
#' @export triSites
#' @author Kristian K Ullrich
triSites<-function(dna, x.pos=NULL, min.ind=0, wlen=25000, start.by=1,
  end.by=NULL, threads=1, pB=FALSE){
    IUPAC_CODE_MAP_LIST<-list(c("A"), c("C"), c("G"), c("T"),
      c("A", "C"), c("A", "G"), c("A", "T"), c("C", "G"), c("C", "T"),
      c("G","T"), c("A", "C", "G"), c("A", "C", "T"), c("A", "G", "T"),
      c("C", "G", "T"), c(), c(), c(), c())
    names(IUPAC_CODE_MAP_LIST)<-c("A", "C", "G", "T", 
      "M", "R", "W", "S", "Y", "K",
      "V", "H", "D", "B",
      "N", "-", "+", ".")
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
        pb<-txtProgressBar(min=0, max=dim(tmp.sw)[2], initial=0, style=3)
    }
    j<-NULL
    registerDoMC(threads)
    OUT<-foreach(j=seq(from=1, to=ncol(tmp.sw)), .combine=rbind) %dopar% {
      START<-tmp.sw[1, j][[1]]
      END<-tmp.sw[2, j][[1]]
      triPOSall<-NA
      triPOS<-NA
      triPOSmasked<-NA
      minIndPOS<-NA
      tmp.seq<-subseq(dna_, START, END)
      cM<-consensusMatrix(tmp.seq)
      minIndPOS<-START - 1 + which(!apply(cM, 2, function(x) {
            sum(x[1:14])
        })>=min.ind)
      if(unique(width(tmp.seq))==1){
          tmp.seq.cM<-t(as.matrix(apply(cM, 1, function(x) {
                ifelse(x>0, 1, 0)
            })))
      }
      if(unique(width(tmp.seq))!=1){
          tmp.seq.cM<-apply(cM, 1, function(x) {
                ifelse(x>0,1,0)
            })
      }
      triPOSall<-START -1 + which(apply(tmp.seq.cM, 1, function(x) {
            length(unique(unlist(unique(
              IUPAC_CODE_MAP_LIST[names(x[x==1])]
            ))))
        })==3)
      triPOS<-triPOSall[which(!triPOSall%in%minIndPOS)]
      triPOSmasked<-triPOSall[which(triPOSall%in%minIndPOS)]
      if(pB){
          setTxtProgressBar(pb, j)
      }
      list(triPOS, triPOSmasked, minIndPOS)
    }
    if(pB){
        setTxtProgressBar(pb, ncol(tmp.sw))
        close(pb)  
    }
    if(is.matrix(OUT)){
        OUT<-list(unlist(OUT[, 1]), unlist(OUT[, 2]), unlist(OUT[, 3]))  
    }
    names(OUT)<-c("triPOS", "triPOSmasked", "minIndPOS")
    names(OUT$triPOS)<-NULL
    names(OUT$triPOSmasked)<-NULL
    names(OUT$minIndPOS)<-NULL
    return(OUT)
}
