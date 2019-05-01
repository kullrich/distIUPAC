#' @title dist.xyStats
#' @name dist.xyStats
#' @description This function calculates \code{distIUPAC} based distances comparing two populations 
#' (x: receiver; y: donor).
#' @importFrom stats as.dist sd
#' @param dIUPAC \code{distIUPAC} distance matrix including \code{sitesUsed} matrix
#' @param x.pos population X positions
#' @param y.pos population Y positions
#' @param x.name population X name
#' @param y.name population Y name
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' #Here, one needs to consider the changed x and y positions due to sub-sampling.
#' #Use the 'xyStats' function which handles this issue automatically
#' CAS.pos_<-1:length(CAS.pos)
#' AFG.pos_<-(length(CAS.pos_)+1):length(c(CAS.pos,AFG.pos))
#' #pairwise deletion
#' CAS_AFG.pairwiseDeletion.dist<-distIUPAC(as.character(MySequences[c(CAS.pos,AFG.pos)]))
#' CAS_AFG.pairwiseDeletion.xyStats<-dist.xyStats(CAS_AFG.pairwiseDeletion.dist, 
#' x.pos=CAS.pos_, y.pos=AFG.pos_, x.name="CAS", y.name="AFG")
#' #global deletion
#' CAS_AFG.globalDeletion.dist<-distIUPAC(as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos)])))
#' CAS_AFG.globalDeletion.xyStats<-dist.xyStats(CAS_AFG.globalDeletion.dist,
#' x.pos=CAS.pos_, y.pos=AFG.pos_, x.name="CAS", y.name="AFG")
#' #compare results
#' rbind(CAS_AFG.pairwiseDeletion.xyStats, CAS_AFG.globalDeletion.xyStats)
#' @export dist.xyStats
#' @author Kristian K Ullrich
dist.xyStats<-function( dIUPAC, x.pos, y.pos, x.name="x", y.name="y"){
  options(scipen=22)
  XNAME<-x.name
  YNAME<-y.name
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
  OUT<-list(XNAME,YNAME,
    dMean.x,dSd.x,dSites.x,dMin.x,dMax.x,
    dMean.y,dSd.y,dSites.y,dMin.y,dMax.y,
    dMean.xy,dSd.xy,dSites.xy,dMin.xy,dMax.xy,
    dTotal.xy,dSweighted.xy,Fst.xy,dRelative.xy)
  names(OUT)<-c("XNAME","YNAME",
    "dMean.x","dSd.x","dSites.x","dMin.x","dMax.x",
    "dMean.y","dSd.y","dSites.y","dMin.y","dMax.y",
    "dMean.xy","dSd.xy","dSites.xy","dMin.xy","dMax.xy",
    "dTotal.xy","dSweighted.xy","Fst.xy","dRelative.xy")
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
  OUT$dMean.xy<-mean(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSd.xy<-sd(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dSites.xy<-mean(dIUPAC$sitesUsed[x.pos,y.pos],na.rm=TRUE)
  OUT$dMin.xy<-min(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dMax.xy<-max(dIUPAC$distIUPAC[x.pos,y.pos],na.rm=TRUE)
  OUT$dTotal.xy<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos,y.pos),c(x.pos,y.pos)]),na.rm=TRUE)
  OUT$dSweighted.xy<-(length(x.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.x + (length(y.pos)/(length(c(x.pos,y.pos)))) * OUT$dMean.y
  OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
  OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
  return(OUT)
}
