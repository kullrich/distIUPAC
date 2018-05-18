#' @title xyMultiStats
#' @name xyMultiStats
#' @description This function calculates \code{distIUPAC} based distances comparing all
#' possible pairwise population combinations.
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param list.pos population positions list
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param threads number of parallel threads
#' @param chr.name chromosome name
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' pop.list<-list(CAS.pos, AFG.pos, SPRE.pos)
#' names(pop.list)<-c("AFG", "CAS", "SPRE")
#' #sliding windows based on base-pair length
#' CAS.AFG.SPRE.xyMultiStats<-xyMultiStats(MySequences, list.pos=pop.list, threads=4)
#' CAS.AFG.SPRE.xyMultiStats
#' #sliding windows based on biSites
#' CAS.AFG.SPRE.xyMultiStats<-xyMultiStats(MySequences, list.pos=pop.list,
#' wtype="biSites", wlen=50, threads=4)
#' CAS.AFG.SPRE.xyMultiStats
#' @export xyMultiStats
#' @author Kristian K Ullrich
xyMultiStats<-function(dna,list.pos,wlen=25000,wjump=25000,wtype="bp",dist="IUPAC",threads=1,chr.name="chr"){
  options(scipen=22)
  if(is.null(names(list.pos))){names(list.pos)<-seq(1,length(list.pos))}
  pop.comb<-combn(length(list.pos),2)
  COMBOUT<-vector("list",length=dim(pop.comb)[2])
  names(COMBOUT)<-apply(pop.comb,2,function(x) paste0(x,collapse="_"))
  for(c.idx in 1:(dim(pop.comb)[2])){
    x.pos<-unlist(list.pos[pop.comb[1,c.idx]])
    y.pos<-unlist(list.pos[pop.comb[2,c.idx]])
    x.name<-names(list.pos[pop.comb[1,c.idx]])
    y.name<-names(list.pos[pop.comb[2,c.idx]])
    dna_<-dna[c(x.pos,y.pos)]
    x.pos_<-seq(1,length(x.pos))
    y.pos_<-seq(length(x.pos_)+1,length(x.pos_)+length(y.pos))
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
    j<-NULL
    pb<-txtProgressBar(min=1,max=dim(tmp.sw)[2],initial=1,style=3)
    registerDoMC(threads)
    OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
      XNAME<-x.name
      YNAME<-y.name
      CHRNAME<-chr.name
      START<-NA
      END<-NA
      dMean.x<-NA
      dSd.x<-NA
      dSites.x<-NA
      dMin.x<-NA
      dMax.x<-NA
      dMean.y<-NA
      dSd.y<-NA
      dSites.y<-NA
      dMin.y<-NA
      dMax.y<-NA
      dMean.xy<-NA
      dSd.xy<-NA
      dSites.xy<-NA
      dMin.xy<-NA
      dMax.xy<-NA
      dTotal.xy<-NA
      dSweighted.xy<-NA
      Fst.xy<-NA
      dRelative.xy<-NA
      OUT<-list(XNAME,YNAME,CHRNAME,START,END,dMean.x,dSd.x,dSites.x,dMin.x,dMax.x,dMean.y,dSd.y,dSites.y,dMin.y,dMax.y,dMean.xy,dSd.xy,dSites.xy,dMin.xy,dMax.xy,dTotal.xy,dSweighted.xy,Fst.xy,dRelative.xy)
      names(OUT)<-c("XNAME","YNAME","CHRNAME","START","END","dMean.x","dSd.x","dSites.x","dMin.x","dMax.x","dMean.y","dSd.y","dSites.y","dMin.y","dMax.y","dMean.xy","dSd.xy","dSites.xy","dMin.xy","dMax.xy","dTotal.xy","dSweighted.xy","Fst.xy","dRelative.xy")
      OUT$START<-tmp.sw[1,j][[1]]
      OUT$END<-tmp.sw[2,j][[1]]
      tmp.seq<-subseq(dna_,OUT$START,OUT$END)
      if(dist=="IUPAC"){
        tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
        OUT$dMean.x<-mean(as.dist(tmp.seq.dist$distIUPAC[x.pos_,x.pos_]),na.rm=T)
        OUT$dSd.x<-sd(as.dist(tmp.seq.dist$distIUPAC[x.pos_,x.pos_]),na.rm=T)
        OUT$dSites.x<-mean(as.dist(tmp.seq.dist$sitesUsed[x.pos_,x.pos_]),na.rm=T)
        OUT$dMin.x<-min(as.dist(tmp.seq.dist$distIUPAC[x.pos_,x.pos_]),na.rm=T)
        OUT$dMax.x<-max(as.dist(tmp.seq.dist$distIUPAC[x.pos_,x.pos_]),na.rm=T)
        OUT$dMean.y<-mean(as.dist(tmp.seq.dist$distIUPAC[y.pos_,y.pos_]),na.rm=T)
        OUT$dSd.y<-sd(as.dist(tmp.seq.dist$distIUPAC[y.pos_,y.pos_]),na.rm=T)
        OUT$dSites.y<-mean(as.dist(tmp.seq.dist$sitesUsed[y.pos_,y.pos_]),na.rm=T)
        OUT$dMin.y<-min(as.dist(tmp.seq.dist$distIUPAC[y.pos_,y.pos_]),na.rm=T)
        OUT$dMax.y<-max(as.dist(tmp.seq.dist$distIUPAC[y.pos_,y.pos_]),na.rm=T)
        OUT$dMean.xy<-mean(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
        OUT$dSd.xy<-sd(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
        OUT$dSites.xy<-mean(tmp.seq.dist$sitesUsed[x.pos_,y.pos_],na.rm=T)
        OUT$dMin.xy<-min(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
        OUT$dMax.xy<-max(tmp.seq.dist$distIUPAC[x.pos_,y.pos_],na.rm=T)
        OUT$dTotal.xy<-mean(as.dist(tmp.seq.dist$distIUPAC),na.rm=T)
        OUT$dSweighted.xy<-(length(x.pos_)/(length(c(x.pos_,y.pos_)))) * OUT$dMean.x + (length(y.pos_)/(length(c(x.pos_,y.pos_)))) * OUT$dMean.y
        OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
        OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
      }
      if(dist!="IUPAC"){
        tmp.seq.dist<-dist.dna(as.DNAbin(dnastring2apealg(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
        tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
        OUT$dMean.x<-mean(as.dist(tmp.seq.dist[x.pos_,x.pos_]),na.rm=T)
        OUT$dSd.x<-sd(as.dist(tmp.seq.dist[x.pos_,x.pos_]),na.rm=T)
        OUT$dSites.x<-mean(as.dist(tmp.seq.sites[x.pos_,x.pos_]),na.rm=T)
        OUT$dMin.x<-min(as.dist(tmp.seq.dist[x.pos_,x.pos_]),na.rm=T)
        OUT$dMax.x<-max(as.dist(tmp.seq.dist[x.pos_,x.pos_]),na.rm=T)
        OUT$dMean.y<-mean(as.dist(tmp.seq.dist[y.pos_,y.pos_]),na.rm=T)
        OUT$dSd.y<-sd(as.dist(tmp.seq.dist[y.pos_,y.pos_]),na.rm=T)
        OUT$dSites.y<-mean(as.dist(tmp.seq.sites[y.pos_,y.pos_]),na.rm=T)
        OUT$dMin.y<-min(as.dist(tmp.seq.dist[y.pos_,y.pos_]),na.rm=T)
        OUT$dMax.y<-max(as.dist(tmp.seq.dist[y.pos_,y.pos_]),na.rm=T)
        OUT$dMean.xy<-mean(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
        OUT$dSd.xy<-sd(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
        OUT$dSites.xy<-mean(tmp.seq.sites[x.pos_,y.pos_],na.rm=T)
        OUT$dMin.xy<-min(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
        OUT$dMax.xy<-max(tmp.seq.dist[x.pos_,y.pos_],na.rm=T)
        OUT$dTotal.xy<-mean(as.dist(tmp.seq.dist),na.rm=T)
        OUT$dSweighted.xy<-(length(x.pos_)/(length(c(x.pos_,y.pos_)))) * OUT$dMean.x + (length(y.pos_)/(length(c(x.pos_,y.pos_)))) * OUT$dMean.y
        OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
        OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
      }
      setTxtProgressBar(pb,j)
      OUT
    }
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)
    COMBOUT[c.idx]<-list(OUT)
  }
  return(COMBOUT)
}
