#' @title xyoStats
#' @name xyoStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations (x: receiver; y: donor).
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param y.pos population Y positions
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param dist distance to use
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param y.name population Y name
#' @param chr.name chromosome name
#' @examples
#' @export xyStats
#' @author Kristian K Ullrich
xyoStats<-function(dna,x.pos,y.pos,wlen=25000,wjump=25000,dist="IUPAC",threads=1,x.name="x",y.name="y",chr.name="chr"){
  options(scipen=22)
  dna_<-dna[c(x.pos,y.pos)]
  x.pos_<-seq(1,length(x.pos))
  y.pos_<-seq(length(x.pos_)+1,length(x.pos_)+length(y.pos))
  tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=1,end.by=unique(width(dna)))
  pb<-txtProgressBar(min=1,max=dim(tmp.sw)[2],initial=1,style=3)
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
    XNAME<-x.name
    YNAME<-y.name
    CHRNAME<-chr.name
    START<-NA
    END<-NA
    dMean.xy<-NA
    dSd.xy<-NA
    dSites.xy<-NA
    dMin.xy<-NA
    OUT<-list(XNAME,YNAME,ONAME,CHRNAME,START,END,dMean.xy,dSd.xy,dSites.xy,dMin.xy)
    names(OUT)<-c("XNAME","YNAME","ONAME","CHRNAME","START","END","dMean.xy","dSd.xy","dSites.xy","dMin.xy")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dSd.xy<-sd(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin(dnastring2apealg(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dSd.xy<-sd(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
    }
    setTxtProgressBar(pb,j)
    OUT
  }
  setTxtProgressBar(pb,dim(tmp.sw)[2])
  close(pb)
  return(OUT)
}