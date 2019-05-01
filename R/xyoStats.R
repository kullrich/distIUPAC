#' @title xyoStats
#' @name xyoStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations
#' (x: receiver; y: donor) with an outgroup population (o: outgroup).
#' In the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2, the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P2 y:P3 o:O] and the populations should be defined
#' for deltaMean and deltaMin statistics as follows [x:P2 y:P3 o:P1].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from P2>>P3, the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P3 y:P2 o:O] and the populations should be defined
#' for deltaMean and deltaMin statistics as follows [x:P3 y:P2 o:P1].
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions [P2 population in the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2]
#' @param y.pos population Y positions [P3 population in the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2]
#' @param o.pos population O positions [P1 or O population in the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2]
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param start.by optional start position
#' @param end.by optional end position
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites with missing data in a global or pairwise way (default is to delete in a global way)
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param y.name population Y name
#' @param o.name population O name
#' @param chr.name chromosome name
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' AFG.SPRE.CAS.xyoStats<-xyoStats(MySequences, x.pos = AFG.pos, y.pos = SPRE.pos, o.pos = CAS.pos,
#' threads = 2, x.name = "AFG", y.name = "SPRE", o.name = "CAS")
#' AFG.SPRE.CAS.xyoStats
#' @export xyoStats
#' @author Kristian K Ullrich
xyoStats<-function(dna, x.pos, y.pos, o.pos, wlen=25000, wjump=25000, start.by=NULL, end.by=NULL, wtype="bp", dist="IUPAC", global.deletion=TRUE, threads=1, x.name="x", y.name="y", o.name="o", chr.name="chr", pB=TRUE){
  options(scipen=22)
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  dna_<-dna[c(x.pos,y.pos,o.pos)]
  x.pos_<-seq(1,length(x.pos))
  y.pos_<-seq(length(x.pos_)+1,length(x.pos_)+length(y.pos))
  o.pos_<-seq(length(x.pos_)+length(y.pos_)+1,length(x.pos_)+length(y.pos_)+length(o.pos))
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=start.by,end.by=end.by)
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna_,c(x.pos_,y.pos_),threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$biPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna_,c(x.pos_,y.pos_),threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$triPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  j<-NULL
  if(pB){
    pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)
  }
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
    RND.xyo<-NA
    Gmin.xyo<-NA
    RNDmin.xyo<-NA
    OUT<-list(XNAME,YNAME,ONAME,CHRNAME,START,END,dMean.xy,dSites.xy,dMin.xy,dMean.xo,dSites.xo,dMin.xo,dMean.yo,dSites.yo,dMin.yo,dMean.xyo,dSites.xyo,dMin.xyo,deltaMean.xyo,deltaMin.xyo,RND.xyo,Gmin.xyo,RNDmin.xyo)
    names(OUT)<-c("XNAME","YNAME","ONAME","CHRNAME","START","END","dMean.xy","dSites.xy","dMin.xy","dMean.xo","dSites.xo","dMin.xo","dMean.yo","dSites.yo","dMin.yo","dMean.xyo","dSites.xyo","dMin.xyo","deltaMean.xyo","deltaMin.xyo","RND.xyo","Gmin.xyo","RNDmin.xyo")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(global.deletion){
      tmp.seq<-globalDeletion(tmp.seq)
    }
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMean.xo<-mean(tmp.seq.dist$distIUPAC[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dSites.xo<-mean(tmp.seq.dist$sitesUsed[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dMin.xo<-min(tmp.seq.dist$distIUPAC[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dMean.yo<-mean(tmp.seq.dist$distIUPAC[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dSites.yo<-mean(tmp.seq.dist$sitesUsed[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dMin.yo<-min(tmp.seq.dist$distIUPAC[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMean.xyo<-mean(c(tmp.seq.dist$distIUPAC[x.pos_,c(y.pos_,o.pos_)],tmp.seq.dist$distIUPAC[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$dSites.xyo<-mean(c(tmp.seq.dist$sitesUsed[x.pos_,c(y.pos_,o.pos_)],tmp.seq.dist$sitesUsed[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$dMin.xyo<-min(c(tmp.seq.dist$distIUPAC[x.pos_,c(y.pos_,o.pos_)],tmp.seq.dist$distIUPAC[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$deltaMean.xyo<-OUT$dMean.xy-OUT$dMean.xo
      OUT$deltaMin.xyo<-OUT$dMin.xy-OUT$dMean.xo
      OUT$RND.xyo<-OUT$dMean.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
      OUT$Gmin.xyo<-OUT$dMin.xy/OUT$dMean.xy
      OUT$RNDmin.xyo<-OUT$dMin.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin(DNAMultipleAlignment(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMean.xo<-mean(tmp.seq.dist[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dSites.xo<-mean(tmp.seq.sites[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dMin.xo<-min(tmp.seq.dist[x.pos_,o.pos_],na.rm=TRUE)
      OUT$dMean.yo<-mean(tmp.seq.dist[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dSites.yo<-mean(tmp.seq.sites[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dMin.yo<-min(tmp.seq.dist[y.pos_,o.pos_],na.rm=TRUE)
      OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=TRUE)
      OUT$dMean.xyo<-mean(c(tmp.seq.dist[x.pos_,c(y.pos_,o.pos_)],tmp.seq.dist[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$dSites.xyo<-mean(c(tmp.seq.sites[x.pos_,c(y.pos_,o.pos_)],tmp.seq.sites[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$dMin.xyo<-min(c(tmp.seq.dist[x.pos_,c(y.pos_,o.pos_)],tmp.seq.dist[y.pos_,o.pos_]),na.rm=TRUE)
      OUT$deltaMean.xyo<-OUT$dMean.xy-OUT$dMean.xo
      OUT$deltaMin.xyo<-OUT$dMin.xy-OUT$dMean.xo
      OUT$RND.xyo<-OUT$dMean.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
      OUT$Gmin.xyo<-OUT$dMin.xy/OUT$dMean.xy
      OUT$RNDmin.xyo<-OUT$dMin.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
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
