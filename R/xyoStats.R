#' @title xyoStats
#' @name xyoStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations (x: receiver; y: donor) with an outgroup population (o: outgroup).
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param y.pos population Y positions
#' @param o.pos population O positions
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param y.name population Y name
#' @param o.name population O name
#' @param chr.name chromosome name
#' @examples
#' @export xyoStats
#' @author Kristian K Ullrich
xyoStats<-function(dna,x.pos,y.pos,o.pos,wlen=25000,wjump=25000,wtype="bp",dist="IUPAC",threads=1,x.name="x",y.name="y",o.name="o",chr.name="chr"){
  options(scipen=22)
  dna_<-dna[c(x.pos,y.pos,o.pos)]
  x.pos_<-seq(1,length(x.pos))
  y.pos_<-seq(length(x.pos_)+1,length(x.pos_)+length(y.pos))
  o.pos_<-seq(length(x.pos_)+length(y.pos_)+1,length(x.pos_)+length(y.pos_)+length(o.pos))
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=1,end.by=unique(width(dna)))  
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna_,c(x.pos_,y.pos_),threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=1,end.by=unique(width(dna)))
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna_,c(x.pos_,y.pos_),threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=1,end.by=unique(width(dna)))
  }
  pb<-txtProgressBar(min=1,max=dim(tmp.sw)[2],initial=1,style=3)
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
    XNAME<-x.name
    YNAME<-y.name
    ONAME<-o.name
    CHRNAME<-chr.name
    START<-NA
    END<-NA
    dMean.xy<-NA
    dSites.xy<-NA
    dMin.xy<-NA
    dMean.xo<-NA
    dSites.xo<-NA
    dMin.xo<-NA
    dMean.yo<-NA
    dSites.yo<-NA
    dMin.yo<-NA
    dMean.xyo<-NA
    dSites.xyo<-NA
    dMin.xyo<-NA
    deltaMean.xyo<-NA
    deltaMin.xyo<-NA
    RNDmin.xyo<-NA
    OUT<-list(XNAME,YNAME,ONAME,CHRNAME,START,END,dMean.xy,dSites.xy,dMin.xy,dMean.xo,dSites.xo,dMin.xo,dMean.yo,dSites.yo,dMin.yo,dMean.xyo,dSites.xyo,dMin.xyo,deltaMean.xyo,deltaMin.xyo,RNDmin.xyo)
    names(OUT)<-c("XNAME","YNAME","ONAME","CHRNAME","START","END","dMean.xy","dSites.xy","dMin.xy","dMean.xo","dSites.xo","dMin.xo","dMean.yo","dSites.yo","dMin.yo","dMean.xyo","dSites.xyo","dMin.xyo","deltaMean.xyo","deltaMin.xyo","RNDmin.xyo")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dMean.xo<-mean(tmp.seq.dist$distIUPAC[x.pos_,o.pos_],na.rm=T)
      OUT$dSites.xo<-mean(tmp.seq.dist$sitesUsed[x.pos_,o.pos_],na.rm=T)
      OUT$dMin.xo<-min(tmp.seq.dist$distIUPAC[x.pos_,o.pos_],na.rm=T)
      OUT$dMean.yo<-mean(tmp.seq.dist$distIUPAC[y.pos_,o.pos_],na.rm=T)
      OUT$dSites.yo<-mean(tmp.seq.dist$sitesUsed[y.pos_,o.pos_],na.rm=T)
      OUT$dMin.yo<-min(tmp.seq.dist$distIUPAC[y.pos_,o.pos_],na.rm=T)
      OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
      OUT$dMean.xyo<-mean(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
      OUT$dSites.xyo<-mean(as.dist(tmp.seq.dist$sitesUsed),na.rm=T)
      OUT$dMin.xyo<-min(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
      OUT$deltaMean.xyo<-OUT$dMean.xy-OUT$dMean.xo
      OUT$deltaMin.xyo<-OUT$dMin.xy-OUT$dMean.xo
      OUT$RNDmin.xyo<-OUT$dMin.xy-((OUT$dMean.xo+OUT$dMean.yo)/2)
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin(dnastring2apealg(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dMean.xo<-mean(tmp.seq.dist[x.pos_,o.pos_],na.rm=T)
      OUT$dSites.xo<-mean(tmp.seq.sites[x.pos_,o.pos_],na.rm=T)
      OUT$dMin.xo<-min(tmp.seq.dist[x.pos_,o.pos_],na.rm=T)
      OUT$dMean.yo<-mean(tmp.seq.dist[y.pos_,o.pos_],na.rm=T)
      OUT$dSites.yo<-mean(tmp.seq.sites[y.pos_,o.pos_],na.rm=T)
      OUT$dMin.yo<-min(tmp.seq.dist[y.pos_,o.pos_],na.rm=T)
      OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=T)
      OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
      OUT$dMean.xyo<-mean(as.dist(tmp.seq.dist),na.rm=T)
      OUT$dSites.xyo<-mean(as.dist(tmp.seq.sites),na.rm=T)
      OUT$dMin.xyo<-min(as.dist(tmp.seq.dist),na.rm=T)
      OUT$deltaMean.xyo<-OUT$dMean.xy-OUT$dMean.xo
      OUT$deltaMin.xyo<-OUT$dMin.xy-OUT$dMean.xo
      OUT$RNDmin.xyo<-OUT$dMin.xy-((OUT$dMean.xo+OUT$dMean.yo)/2)
    }
    setTxtProgressBar(pb,j)
    OUT
  }
  setTxtProgressBar(pb,dim(tmp.sw)[2])
  close(pb)
  return(OUT)
}
