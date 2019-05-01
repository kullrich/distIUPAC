#' @title xStats
#' @name xStats
#' @description This function calculates \code{distIUPAC} based distances within one population
#' (x: receiver; x: donor).
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param start.by optional start position
#' @param end.by optional end position
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites with missing data in a global or pairwise way (default is to delete in a global way)
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param chr.name chromosome name
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' CAS.xStats<-xStats(MySequences, x.pos = CAS.pos, x.name = "CAS", threads = 2)
#' CAS.xStats
#' @export xStats
#' @author Kristian K Ullrich
xStats<-function(dna, x.pos=NULL, wlen=25000, wjump=25000, start.by=NULL, end.by=NULL, wtype="bp", dist="IUPAC", global.deletion=TRUE, threads=1, x.name="x", chr.name="chr", pB=TRUE){
  options(scipen=22)
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  if(is.null(x.pos)){
    x.pos<-seq(1,length(dna))
  }
  dna_<-dna[x.pos]
  x.pos_<-seq(1,length(x.pos))
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=start.by,end.by=end.by)
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$biPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$triPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  j<-NULL
  if(pB){
    pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)
  }
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
    XNAME<-x.name
    CHRNAME<-chr.name
    START<-NA
    END<-NA
    dMean.x<-NA
    dSd.x<-NA
    dSites.x<-NA
    dMin.x<-NA
    dMax.x<-NA
    OUT<-list(XNAME,CHRNAME,START,END,dMean.x,dSd.x,dSites.x,dMin.x,dMax.x)
    names(OUT)<-c("XNAME","CHRNAME","START","END","dMean.x","dSd.x","dSites.x","dMin.x","dMax.x")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(global.deletion){
      tmp.seq<-globalDeletion(tmp.seq)
    }
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dMean.x<-mean(as.dist(tmp.seq.dist$distIUPAC),na.rm=TRUE)
      OUT$dSd.x<-sd(as.dist(tmp.seq.dist$distIUPAC),na.rm=TRUE)
      OUT$dSites.x<-mean(as.dist(tmp.seq.dist$sitesUsed),na.rm=TRUE)
      OUT$dMin.x<-min(as.dist(tmp.seq.dist$distIUPAC),na.rm=TRUE)
      OUT$dMax.x<-max(as.dist(tmp.seq.dist$distIUPAC),na.rm=TRUE)
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin(DNAMultipleAlignment(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dMean.x<-mean(as.dist(tmp.seq.dist),na.rm=TRUE)
      OUT$dSd.x<-sd(as.dist(tmp.seq.dist),na.rm=TRUE)
      OUT$dSites.x<-mean(as.dist(tmp.seq.sites),na.rm=TRUE)
      OUT$dMin.x<-min(as.dist(tmp.seq.dist),na.rm=TRUE)
      OUT$dMax.x<-max(as.dist(tmp.seq.dist),na.rm=TRUE)
    }
    if(pB){
      setTxtProgressBar(pb,j)
    }
    OUT
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)
  }
  return(OUT)
}
distIUPAC2xStats<-function( seq.distIUPAC, x.name="x", x.pos=NULL){
  options(scipen=22)
  if(is.null(x.pos)){
    x.pos<-1:dim(seq.distIUPAC$distIUPAC)[1]
  }
  XNAME<-x.name
  dMean.x<-NA
  dSd.x<-NA
  dSites.x<-NA
  dMin.x<-NA
  dMax.x<-NA
  OUT<-list(XNAME,dMean.x,dSd.x,dSites.x,dMin.x,dMax.x)
  names(OUT)<-c("XNAME","dMean.x","dSd.x","dSites.x","dMin.x","dMax.x")
  OUT$dMean.x<-mean(as.dist(seq.distIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSd.x<-sd(as.dist(seq.distIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSites.x<-mean(as.dist(seq.distIUPAC$sitesUsed[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMin.x<-min(as.dist(seq.distIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMax.x<-max(as.dist(seq.distIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  return(OUT)
}
