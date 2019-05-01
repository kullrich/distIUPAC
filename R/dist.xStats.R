#' @title dist.xStats
#' @name dist.xStats
#' @description This function calculates \code{distIUPAC} based distances within one population
#' (x: receiver; x: donor).
#' @importFrom stats as.dist sd
#' @param dIUPAC \code{distIUPAC} distance matrix including \code{sitesUsed} matrix
#' @param x.pos population X positions
#' @param x.name population X name
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' #pairwise deletion
#' CAS.pairwiseDeletion.dist<-distIUPAC(as.character(MySequences[CAS.pos]))
#' CAS.pairwiseDeletion.xStats<-dist.xStats(CAS.pairwiseDeletion.dist, x.name="CAS")
#' #global deletion
#' CAS.globalDeletion.dist<-distIUPAC(as.character(globalDeletion(MySequences[CAS.pos])))
#' CAS.globalDeletion.xStats<-dist.xStats(CAS.globalDeletion.dist, x.name="CAS")
#' #compare results
#' rbind(CAS.pairwiseDeletion.xStats, CAS.globalDeletion.xStats)
#' @export dist.xStats
#' @author Kristian K Ullrich
dist.xStats<-function( dIUPAC, x.pos=NULL, x.name="x"){
  options(scipen=22)
  if(is.null(x.pos)){
    x.pos<-1:dim(dIUPAC$distIUPAC)[1]
  }
  XNAME<-x.name
  dMean.x<-NA
  dSd.x<-NA
  dSites.x<-NA
  dMin.x<-NA
  dMax.x<-NA
  OUT<-list(XNAME,dMean.x,dSd.x,dSites.x,dMin.x,dMax.x)
  names(OUT)<-c("XNAME","dMean.x","dSd.x","dSites.x","dMin.x","dMax.x")
  OUT$dMean.x<-mean(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSd.x<-sd(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dSites.x<-mean(as.dist(dIUPAC$sitesUsed[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMin.x<-min(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  OUT$dMax.x<-max(as.dist(dIUPAC$distIUPAC[x.pos,x.pos]),na.rm=TRUE)
  return(OUT)
}
