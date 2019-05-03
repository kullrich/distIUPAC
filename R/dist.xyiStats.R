#' @title dist.xyiStats
#' @name dist.xyiStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations
#' (x: receiver; y: donor) with an ingroup population (i: ingroup).
#' In a four-taxin scenario (((P1,P2),P3),O) with geneflow from P3>>P2, the populations should be defined
#' for deltaMean and deltaMin statistics as follows [x:P2 y:P3 i:P1].
#' Accordingly in the four-taxin scenario (((P1,P2),P3),O) with geneflow from P2>>P3, the populations should be defined
#' for deltaMean and deltaMin statistics as follows [x:P3 y:P2 i:P1].
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
#' IRA.pos<-71:78
#' #Here, one needs to consider the changed x, y and i positions due to sub-sampling.
#' #Use the 'xyiStats' function which handles this issue automatically
#' CAS.pos_<-1:length(CAS.pos)
#' AFG.pos_<-(length(CAS.pos_)+1):length(c(CAS.pos,AFG.pos))
#' IRA.pos_<-(length(c(CAS.pos,AFG.pos))+1):length(c(CAS.pos,AFG.pos,IRA.pos))
#' #pairwise deletion
#' CAS_AFG_IRA.pairwiseDeletion.dist<-distIUPAC(as.character(MySequences[c(CAS.pos,AFG.pos,IRA.pos)]))
#' CAS_AFG_IRA.pairwiseDeletion.xyiStats<-dist.xyiStats(CAS_AFG_IRA.pairwiseDeletion.dist, 
#' x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_, x.name="CAS", y.name="AFG", i.name="IRA")
#' #global deletion
#' CAS_AFG_IRA.globalDeletion.dist<-distIUPAC(as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos,IRA.pos)])))
#' CAS_AFG_IRA.globalDeletion.xyiStats<-dist.xyiStats(CAS_AFG_IRA.globalDeletion.dist,
#' x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_, x.name="CAS", y.name="AFG", i.name="IRA")
#' #compare results
#' rbind(CAS_AFG_IRA.pairwiseDeletion.xyiStats, CAS_AFG_IRA.globalDeletion.xyiStats)
#' @export dist.xyiStats
#' @author Kristian K Ullrich
dist.xyiStats<-function( dIUPAC, x.pos, y.pos, i.pos, x.name="x", y.name="y", i.name="i"){
  options(scipen=22)
  XNAME<-x.name
  YNAME<-y.name
  INAME<-i.name
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
  dMean.i<-NA
  dSd.i<-NA
  dSites.i<-NA
  dMin.i<-NA
  dMax.i<-NA
  dMean.xy<-NA
  dSd.xy<-NA
  dSites.xy<-NA
  dMin.xy<-NA
  dMax.xy<-NA
  dTotal.xy<-NA
  dSweighted.xy<-NA
  Fst.xy<-NA
  dRelative.xy<-NA
  dMean.xi<-NA
  dSd.xi<-NA
  dSites.xi<-NA
  dMin.xi<-NA
  dMax.xi<-NA
  dTotal.xi<-NA
  dSweighted.xi<-NA
  Fst.xi<-NA
  dRelative.xi<-NA
  dMean.yi<-NA
  dSd.yi<-NA
  dSites.yi<-NA
  dMin.yi<-NA
  dMax.yi<-NA
  dTotal.yi<-NA
  dSweighted.yi<-NA
  Fst.yi<-NA
  dRelative.yi<-NA
  dMean.xyi<-NA
  dSd.xyi<-NA
  dSites.xyi<-NA
  dMin.xyi<-NA
  dMax.xyi<-NA
  dTotal.xyi<-NA
  deltaMean.xyi<-NA
  deltaMin.xyi<-NA
  RND.xyi<-NA
  Gmin.xyi<-NA
  RNDmin.xyi<-NA
  OUT<-list(XNAME,YNAME,INAME,
    dMean.x,dSd.x,dSites.x,dMin.x,dMax.x,
    dMean.y,dSd.y,dSites.y,dMin.y,dMax.y,
    dMean.i,dSd.i,dSites.i,dMin.i,dMax.i,
    dMean.xy,dSd.xy,dSites.xy,dMin.xy,dMax.xy,
    dTotal.xy,dSweighted.xy,Fst.xy,dRelative.xy,
    dMean.xi,dSd.xi,dSites.xi,dMin.xi,dMax.xi,
    dTotal.xi,dSweighted.xi,Fst.xi,dRelative.xi,
    dMean.yi,dSd.yi,dSites.yi,dMin.yi,dMax.yi,
    dTotal.yi,dSweighted.yi,Fst.yi,dRelative.yi,
    dMean.xyi,dSd.xyi,dSites.xyi,dMin.xyi,dMax.xyi,dTotal.xyi,
    RND.xyi,Gmin.xyi,RNDmin.xyi)
  names(OUT)<-c("XNAME","YNAME","ONAME",
    "dMean.x","dSd.x","dSites.x","dMin.x","dMax.x",
    "dMean.y","dSd.y","dSites.y","dMin.y","dMax.y",
    "dMean.i","dSd.i","dSites.i","dMin.i","dMax.i",
    "dMean.xy","dSd.xy","dSites.xy","dMin.xy","dMax.xy",
    "dTotal.xy","dSweighted.xy","Fst.xy","dRelative.xy",
    "dMean.xi","dSd.xi","dSites.xi","dMin.xi","dMax.xi",
    "dTotal.xi","dSweighted.xi","Fst.xi","dRelative.xi",
    "dMean.yi","dSd.yi","dSites.yi","dMin.yi","dMax.yi",
    "dTotal.yi","dSweighted.yi","Fst.yi","dRelative.yi",
    "dMean.xyi","dSd.xyi","dSites.xyi","dMin.xyi","dMax.xyi","dTotal.xyi",
    "RND.xyi","Gmin.xyi","RNDmin.xyi")
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
  OUT$dMean.i<-mean(as.dist(dIUPAC$distIUPAC[i.pos,i.pos]),na.rm=TRUE)
  OUT$dSd.i<-sd(as.dist(dIUPAC$distIUPAC[i.pos,i.pos]),na.rm=TRUE)
  OUT$dSites.i<-mean(as.dist(dIUPAC$sitesUsed[i.pos,i.pos]),na.rm=TRUE)
  OUT$dMin.i<-min(as.dist(dIUPAC$distIUPAC[i.pos,i.pos]),na.rm=TRUE)
  OUT$dMax.i<-max(as.dist(dIUPAC$distIUPAC[i.pos,i.pos]),na.rm=TRUE)
  OUT$dMean.xy<-mean(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSd.xy<-sd(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSites.xy<-mean(dIUPAC$sitesUsed[x.pos,y.pos],na.rm=TRUE)
  OUT$dMin.xy<-min(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dMax.xy<-max(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dTotal.xy<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,y.pos),c(x.pos,y.pos)]),na.rm=TRUE)
  OUT$dSweighted.xy<-(length(x.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.x + (length(y.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.y
  OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
  OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
  OUT$dMean.xi<-mean(dIUPAC$distIUPAC[x.pos,i.pos],na.rm=TRUE)
  OUT$dSd.xi<-sd(dIUPAC$distIUPAC[x.pos,i.pos],na.rm=TRUE)
  OUT$dSites.xi<-mean(dIUPAC$sitesUsed[x.pos,i.pos],na.rm=TRUE)
  OUT$dMin.xi<-min(dIUPAC$distIUPAC[x.pos,i.pos],na.rm=TRUE)
  OUT$dMax.xi<-max(dIUPAC$distIUPAC[x.pos,i.pos],na.rm=TRUE)
  OUT$dTotal.xi<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,i.pos),c(x.pos,i.pos)]),na.rm=TRUE)
  OUT$dSweighted.xi<-(length(x.pos)/(length(c(x.pos,i.pos)))) * OUT$dMean.x + (length(i.pos)/(length(c(x.pos,i.pos)))) * OUT$dMean.i
  OUT$Fst.xi<-(OUT$dTotal.xi- OUT$dSweighted.xi) / OUT$dTotal.xi
  OUT$dRelative.xi<-OUT$dMean.xi - OUT$dSweighted.xi
  OUT$dMean.yi<-mean(dIUPAC$distIUPAC[y.pos,i.pos],na.rm=TRUE)
  OUT$dSd.yi<-sd(dIUPAC$distIUPAC[y.pos,i.pos],na.rm=TRUE)
  OUT$dSites.yi<-mean(dIUPAC$sitesUsed[y.pos,i.pos],na.rm=TRUE)
  OUT$dMin.yi<-min(dIUPAC$distIUPAC[y.pos,i.pos],na.rm=TRUE)
  OUT$dMax.yi<-max(dIUPAC$distIUPAC[y.pos,i.pos],na.rm=TRUE)
  OUT$dTotal.yi<-mean(as.dist(dIUPAC$distIUPAC[c(y.pos,i.pos),c(y.pos,i.pos)]),na.rm=TRUE)
  OUT$dSweighted.yi<-(length(y.pos)/(length(c(y.pos,i.pos)))) * OUT$dMean.y + (length(i.pos)/(length(c(x.pos,i.pos)))) * OUT$dMean.i
  OUT$Fst.yi<-(OUT$dTotal.yi- OUT$dSweighted.yi) / OUT$dTotal.yi
  OUT$dRelative.yi<-OUT$dMean.yi - OUT$dSweighted.yi
  OUT$dMean.xyi<-mean(c(dIUPAC$distIUPAC[x.pos,c(y.pos,i.pos)],dIUPAC$distIUPAC[y.pos,i.pos]),na.rm=TRUE)
  OUT$dSd.xyi<-sd(c(dIUPAC$distIUPAC[x.pos,c(y.pos,i.pos)],dIUPAC$distIUPAC[y.pos,i.pos]),na.rm=TRUE)
  OUT$dSites.xyi<-mean(c(dIUPAC$sitesUsed[x.pos,c(y.pos,i.pos)],dIUPAC$sitesUsed[y.pos,i.pos]),na.rm=TRUE)
  OUT$dMin.xyi<-min(c(dIUPAC$distIUPAC[x.pos,c(y.pos,i.pos)],dIUPAC$distIUPAC[y.pos,i.pos]),na.rm=TRUE)
  OUT$dMax.xyi<-max(c(dIUPAC$distIUPAC[x.pos,c(y.pos,i.pos)],dIUPAC$distIUPAC[y.pos,i.pos]),na.rm=TRUE)
  OUT$dTotal.xyi<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,y.pos,i.pos),c(x.pos,y.pos,i.pos)]),na.rm=TRUE)
  OUT$deltaMean.xyi<-OUT$dMean.xy-OUT$dMean.xi
  OUT$deltaMin.xyi<-OUT$dMin.xy-OUT$dMean.xi
  OUT$RND.xyi<-OUT$dMean.xy/((OUT$dMean.xi+OUT$dMean.yi)/2)
  OUT$Gmin.xyi<-OUT$dMin.xy/OUT$dMean.xy
  OUT$RNDmin.xyi<-OUT$dMin.xy/((OUT$dMean.xi+OUT$dMean.yi)/2)
  return(OUT)
}
