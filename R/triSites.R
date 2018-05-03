#' @title triSites
#' @name triSites
#' @description This function returns tri-allelic site positions given a dna object \code{DNAStringSet}, also with IUPAC code.
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param wlen sliding windows length
#' @param threads number of parallel threads
#' @examples
#' @export triSites
#' @author Kristian K Ullrich
IUPAC_CODE_MAP_LIST<-list(c("A"),c("C"),c("G"),c("T"),c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"),c("A","C","G"),c("A","C","T"),c("A","G","T"),c("C","G","T"),c(),c(),c(),c())
names(IUPAC_CODE_MAP_LIST)<-c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B","N","-","+",".")
triSites<-function(dna,x.pos,wlen=25000,threads=1,pB=TRUE){
  options(scipen=22)
  dna_<-dna[x.pos]
  tmp.sw<-swgen(wlen=wlen,wjump=wlen,start.by=1,end.by=unique(width(dna)))
  if(pB){
    pb<-txtProgressBar(min=1,max=dim(tmp.sw)[2],initial=1,style=3)
  }
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=c) %dopar% {
    START<-NA
    END<-NA
    triPOS<-NA
    OUT<-list(START,END,triPOS)
    names(OUT)<-c("START","END","triPOS")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    tmp.seq.cM<-apply(consensusMatrix(tmp.seq),1,function(x) ifelse(x>0,1,0))
    triPOS<-OUT$START-1+which(apply(tmp.seq.cM,1,function(x) length(unlist(unique(IUPAC_CODE_MAP_LIST[names(x[x==1])]))))==3)
    if(pB){
      setTxtProgressBar(pb,j)
    }
    triPOS
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)  
  }
  return(OUT)
}