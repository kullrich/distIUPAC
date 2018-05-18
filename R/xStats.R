#' @title xStats
#' @name xStats
#' @description This function calculates \code{distIUPAC} based distances within one population
#' (x: receiver; x: donor).
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param chr.name chromosome name
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' CAS.xStats<-xStats(MySequences, x.pos=CAS.pos, x.name="CAS")
#' CAS.xStats
#' @export xStats
#' @author Kristian K Ullrich
xStats<-function(dna, x.pos=NULL, wlen=25000, wjump=25000, wtype="bp", dist="IUPAC", threads=1, x.name="x", chr.name="chr"){
  options(scipen=22)
  if(is.null(x.pos)){dna_<-dna}
  if(!is.null(x.pos)){dna_<-dna[x.pos]}
  x.pos_<-seq(1,length(x.pos))
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=1,end.by=unique(width(dna)))  
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=1,end.by=unique(width(dna)))
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=1,end.by=unique(width(dna)))
  }
  j<-NULL
  pb<-txtProgressBar(min=1,max=dim(tmp.sw)[2],initial=1,style=3)
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
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dMean.x<-mean(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
      OUT$dSd.x<-sd(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
      OUT$dSites.x<-mean(as.dist(tmp.seq.dist$sitesUsed),na.rm=T)
      OUT$dMin.x<-min(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
      OUT$dMax.x<-max(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin(dnastring2apealg(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dMean.x<-mean(as.dist(tmp.seq.dist),na.rm=T)
      OUT$dSd.x<-sd(as.dist(tmp.seq.dist),na.rm=T)
      OUT$dSites.x<-mean(as.dist(tmp.seq.sites),na.rm=T)
      OUT$dMin.x<-min(as.dist(tmp.seq.dist),na.rm=T)
      OUT$dMax.x<-max(as.dist(tmp.seq.dist),na.rm=T)
    }
    setTxtProgressBar(pb,j)
    OUT
  }
  setTxtProgressBar(pb,dim(tmp.sw)[2])
  close(pb)
  return(OUT)
}
