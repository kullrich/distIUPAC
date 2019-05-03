#' @title dist.xyiStats
#' @name dist.xyiStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations
#' (x: receiver; y: donor) with an ingroup population (i: ingroup).
#' In a four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2, the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P2 y:P3 i:O].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from P2>>P3, the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P3 y:P2 i:O].
#' @importFrom stats as.dist sd
#' @param dIUPAC \code{distIUPAC} distance matrix including \code{sitesUsed} matrix
#' @param x.pos population X positions
#' @param y.pos population Y positions
#' @param i.pos population I positions
#' @param x.name population X name
#' @param y.name population Y name
#' @param i.name population I name
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' #Here, one needs to consider the changed x, y and o positions due to sub-sampling.
#' #Use the 'xyoStats' function which handles this issue automatically
#' CAS.pos_<-1:length(CAS.pos)
#' AFG.pos_<-(length(CAS.pos_)+1):length(c(CAS.pos,AFG.pos))
#' SPRE.pos_<-(length(c(CAS.pos,AFG.pos))+1):length(c(CAS.pos,AFG.pos,SPRE.pos))
#' #pairwise deletion
#' CAS_AFG_SPRE.pairwiseDeletion.dist<-distIUPAC(as.character(MySequences[c(CAS.pos,AFG.pos,SPRE.pos)]))
#' CAS_AFG_SPRE.pairwiseDeletion.xyoStats<-dist.xyoStats(CAS_AFG_SPRE.pairwiseDeletion.dist, 
#' x.pos=CAS.pos_, y.pos=AFG.pos_, o.pos=SPRE.pos_, x.name="CAS", y.name="AFG", o.name="SPRE")
#' #global deletion
#' CAS_AFG_SPRE.globalDeletion.dist<-distIUPAC(as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos,SPRE.pos)])))
#' CAS_AFG_SPRE.globalDeletion.xyoStats<-dist.xyoStats(CAS_AFG_SPRE.globalDeletion.dist,
#' x.pos=CAS.pos_, y.pos=AFG.pos_, o.pos=SPRE.pos_, x.name="CAS", y.name="AFG", o.name="SPRE")
#' #compare results
#' rbind(CAS_AFG_SPRE.pairwiseDeletion.xyoStats, CAS_AFG_SPRE.globalDeletion.xyoStats)
#' @export dist.xyoStats
#' @author Kristian K Ullrich
dist.xyoStats<-function( dIUPAC, x.pos, y.pos, o.pos, x.name="x", y.name="y", o.name="o"){
  options(scipen=22)
  XNAME<-x.name
  YNAME<-y.name
  ONAME<-o.name
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
  dMean.o<-NA
  dSd.o<-NA
  dSites.o<-NA
  dMin.o<-NA
  dMax.o<-NA
  dMean.xy<-NA
  dSd.xy<-NA
  dSites.xy<-NA
  dMin.xy<-NA
  dMax.xy<-NA
  dTotal.xy<-NA
  dSweighted.xy<-NA
  Fst.xy<-NA
  dRelative.xy<-NA
  dMean.xo<-NA
  dSd.xo<-NA
  dSites.xo<-NA
  dMin.xo<-NA
  dMax.xo<-NA
  dTotal.xo<-NA
  dSweighted.xo<-NA
  Fst.xo<-NA
  dRelative.xo<-NA
  dMean.yo<-NA
  dSd.yo<-NA
  dSites.yo<-NA
  dMin.yo<-NA
  dMax.yo<-NA
  dTotal.yo<-NA
  dSweighted.yo<-NA
  Fst.yo<-NA
  dRelative.yo<-NA
  dMean.xyo<-NA
  dSd.xyo<-NA
  dSites.xyo<-NA
  dMin.xyo<-NA
  dMax.xyo<-NA
  dTotal.xyo<-NA
  deltaMean.xyo<-NA
  deltaMin.xyo<-NA
  RND.xyo<-NA
  Gmin.xyo<-NA
  RNDmin.xyo<-NA
  OUT<-list(XNAME,YNAME,ONAME,
    dMean.x,dSd.x,dSites.x,dMin.x,dMax.x,
    dMean.y,dSd.y,dSites.y,dMin.y,dMax.y,
    dMean.o,dSd.o,dSites.o,dMin.o,dMax.o,
    dMean.xy,dSd.xy,dSites.xy,dMin.xy,dMax.xy,
    dTotal.xy,dSweighted.xy,Fst.xy,dRelative.xy,
    dMean.xo,dSd.xo,dSites.xo,dMin.xo,dMax.xo,
    dTotal.xo,dSweighted.xo,Fst.xo,dRelative.xo,
    dMean.yo,dSd.yo,dSites.yo,dMin.yo,dMax.yo,
    dTotal.yo,dSweighted.yo,Fst.yo,dRelative.yo,
    dMean.xyo,dSd.xyo,dSites.xyo,dMin.xyo,dMax.xyo,dTotal.xyo,
    deltaMean.xyo,deltaMin.xyo,
    RND.xyo,Gmin.xyo,RNDmin.xyo)
  names(OUT)<-c("XNAME","YNAME","ONAME",
    "dMean.x","dSd.x","dSites.x","dMin.x","dMax.x",
    "dMean.y","dSd.y","dSites.y","dMin.y","dMax.y",
    "dMean.o","dSd.o","dSites.o","dMin.o","dMax.o",
    "dMean.xy","dSd.xy","dSites.xy","dMin.xy","dMax.xy",
    "dTotal.xy","dSweighted.xy","Fst.xy","dRelative.xy",
    "dMean.xo","dSd.xo","dSites.xo","dMin.xo","dMax.xo",
    "dTotal.xo","dSweighted.xo","Fst.xo","dRelative.xo",
    "dMean.yo","dSd.yo","dSites.yo","dMin.yo","dMax.yo",
    "dTotal.yo","dSweighted.yo","Fst.yo","dRelative.yo",
    "dMean.xyo","dSd.xyo","dSites.xyo","dMin.xyo","dMax.xyo","dTotal.xyo",
    "deltaMean.xyo","deltaMin.xyo",
    "RND.xyo","Gmin.xyo","RNDmin.xyo")
  OUT$dMean.x<-mean(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSd.x<-sd(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSites.x<-mean(as.dist(dIUPAC$sitesUsed[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMin.x<-min(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMax.x<-max(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMean.y<-mean(as.dist(dIUPAC$distIUPAC[y.pos,y.pos]),na.rm=TRUE)
  OUT$dSd.y<-sd(as.dist(dIUPAC$distIUPAC[y.pos,y.pos]),na.rm=TRUE)
  OUT$dSites.y<-mean(as.dist(dIUPAC$sitesUsed[y.pos,y.pos]),na.rm=TRUE)
  OUT$dMin.y<-min(as.dist(dIUPAC$distIUPAC[y.pos,y.pos]),na.rm=TRUE)
  OUT$dMax.y<-max(as.dist(dIUPAC$distIUPAC[y.pos,y.pos]),na.rm=TRUE)
  OUT$dMean.o<-mean(as.dist(dIUPAC$distIUPAC[o.pos,o.pos]),na.rm=TRUE)
  OUT$dSd.o<-sd(as.dist(dIUPAC$distIUPAC[o.pos,o.pos]),na.rm=TRUE)
  OUT$dSites.o<-mean(as.dist(dIUPAC$sitesUsed[o.pos,o.pos]),na.rm=TRUE)
  OUT$dMin.o<-min(as.dist(dIUPAC$distIUPAC[o.pos,o.pos]),na.rm=TRUE)
  OUT$dMax.o<-max(as.dist(dIUPAC$distIUPAC[o.pos,o.pos]),na.rm=TRUE)
  OUT$dMean.xy<-mean(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSd.xy<-sd(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSites.xy<-mean(dIUPAC$sitesUsed[x.pos,y.pos],na.rm=TRUE)
  OUT$dMin.xy<-min(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dMax.xy<-max(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dTotal.xy<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,y.pos),c(x.pos,y.pos)]),na.rm=TRUE)
  OUT$dSweighted.xy<-(length(x.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.x + (length(y.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.y
  OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
  OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
  OUT$dMean.xo<-mean(dIUPAC$distIUPAC[x.pos,o.pos],na.rm=TRUE)
  OUT$dSd.xo<-sd(dIUPAC$distIUPAC[x.pos,o.pos],na.rm=TRUE)
  OUT$dSites.xo<-mean(dIUPAC$sitesUsed[x.pos,o.pos],na.rm=TRUE)
  OUT$dMin.xo<-min(dIUPAC$distIUPAC[x.pos,o.pos],na.rm=TRUE)
  OUT$dMax.xo<-max(dIUPAC$distIUPAC[x.pos,o.pos],na.rm=TRUE)
  OUT$dTotal.xo<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,o.pos),c(x.pos,o.pos)]),na.rm=TRUE)
  OUT$dSweighted.xo<-(length(x.pos)/(length(c(x.pos,o.pos)))) * OUT$dMean.x + (length(o.pos)/(length(c(x.pos,o.pos)))) * OUT$dMean.o
  OUT$Fst.xo<-(OUT$dTotal.xo- OUT$dSweighted.xo) / OUT$dTotal.xo
  OUT$dRelative.xo<-OUT$dMean.xo - OUT$dSweighted.xo
  OUT$dMean.yo<-mean(dIUPAC$distIUPAC[y.pos,o.pos],na.rm=TRUE)
  OUT$dSd.yo<-sd(dIUPAC$distIUPAC[y.pos,o.pos],na.rm=TRUE)
  OUT$dSites.yo<-mean(dIUPAC$sitesUsed[y.pos,o.pos],na.rm=TRUE)
  OUT$dMin.yo<-min(dIUPAC$distIUPAC[y.pos,o.pos],na.rm=TRUE)
  OUT$dMax.yo<-max(dIUPAC$distIUPAC[y.pos,o.pos],na.rm=TRUE)
  OUT$dTotal.yo<-mean(as.dist(dIUPAC$distIUPAC[c(y.pos,o.pos),c(y.pos,o.pos)]),na.rm=TRUE)
  OUT$dSweighted.yo<-(length(y.pos)/(length(c(y.pos,o.pos)))) * OUT$dMean.y + (length(o.pos)/(length(c(x.pos,o.pos)))) * OUT$dMean.o
  OUT$Fst.yo<-(OUT$dTotal.yo- OUT$dSweighted.yo) / OUT$dTotal.yo
  OUT$dRelative.yo<-OUT$dMean.yo - OUT$dSweighted.yo
  OUT$dMean.xyo<-mean(c(dIUPAC$distIUPAC[x.pos,c(y.pos,o.pos)],dIUPAC$distIUPAC[y.pos,o.pos]),na.rm=TRUE)
  OUT$dSd.xyo<-sd(c(dIUPAC$distIUPAC[x.pos,c(y.pos,o.pos)],dIUPAC$distIUPAC[y.pos,o.pos]),na.rm=TRUE)
  OUT$dSites.xyo<-mean(c(dIUPAC$sitesUsed[x.pos,c(y.pos,o.pos)],dIUPAC$sitesUsed[y.pos,o.pos]),na.rm=TRUE)
  OUT$dMin.xyo<-min(c(dIUPAC$distIUPAC[x.pos,c(y.pos,o.pos)],dIUPAC$distIUPAC[y.pos,o.pos]),na.rm=TRUE)
  OUT$dMax.xyo<-max(c(dIUPAC$distIUPAC[x.pos,c(y.pos,o.pos)],dIUPAC$distIUPAC[y.pos,o.pos]),na.rm=TRUE)
  OUT$dTotal.xyo<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,y.pos,o.pos),c(x.pos,y.pos,o.pos)]),na.rm=TRUE)
  OUT$deltaMean.xyo<-OUT$dMean.xy-OUT$dMean.xo
  OUT$deltaMin.xyo<-OUT$dMin.xy-OUT$dMean.xo
  OUT$RND.xyo<-OUT$dMean.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
  OUT$Gmin.xyo<-OUT$dMin.xy/OUT$dMean.xy
  OUT$RNDmin.xyo<-OUT$dMin.xy/((OUT$dMean.xo+OUT$dMean.yo)/2)
  return(OUT)
}
