#' @title dist.xyiStats
#' @name dist.xyiStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver; y: donor) with an ingroup
#' population (i: ingroup).
#' In a four-taxin scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined for deltaMean and deltaMin statistics
#' as follows [x:P2 y:P3 i:P1].
#' Accordingly in the four-taxin scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined for deltaMean and deltaMin
#' statistics as follows [x:P3 y:P2 i:P1].
#' @importFrom stats as.dist sd setNames
#' @param tmpSEQ \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [mandatory]
#' @param y.pos population Y positions [mandatory]
#' @param i.pos population I positions [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param i.name population I name [default: "i"]
#' @param dist distance to use [default: IUPAC] or choose one model as in
#' \link[ape]{dist.dna} [default: "IUPAC"]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC} [default: 1]
#' @seealso \code{\link[distIUPAC]{xyiStats}},
#' \code{\link[distIUPAC]{distIUPAC}}, \code{\link[distIUPAC]{rcpp_distIUPAC}}
#' @examples
#' ##Use the 'xyiStats' function to handle population assignment automatically
#' 
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' IRA.pos<-71:78
#' 
#' ##Here, one needs to consider the changed x, y and i positions due
#' ##to sub-sampling.
#' CAS.pos_<-seq_along(CAS.pos)
#' AFG.pos_<-seq(from=length(CAS.pos_)+1, to=length(c(CAS.pos, AFG.pos)))
#' IRA.pos_<-seq(from=length(c(CAS.pos, AFG.pos))+1,
#' to=length(c(CAS.pos, AFG.pos, IRA.pos)))
#' 
#' ##pairwise deletion
#' CAS_AFG_IRA.pairwiseDeletion.dist<-distIUPAC(
#'   as.character(MySequences[c(CAS.pos,AFG.pos,IRA.pos)]))
#' CAS_AFG_IRA.pairwiseDeletion.xyiStats<-dist.xyiStats(
#'   CAS_AFG_IRA.pairwiseDeletion.dist, 
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_,
#'   x.name="CAS", y.name="AFG", i.name="IRA")
#' 
#' ##global deletion
#' CAS_AFG_IRA.globalDeletion.dist<-distIUPAC(
#'   as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos,IRA.pos)])))
#' CAS_AFG_IRA.globalDeletion.xyiStats<-dist.xyiStats(
#'   CAS_AFG_IRA.globalDeletion.dist,
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_,
#'   x.name="CAS", y.name="AFG", i.name="IRA")
#' 
#' ##compare results
#' rbind(CAS_AFG_IRA.pairwiseDeletion.xyiStats,
#'   CAS_AFG_IRA.globalDeletion.xyiStats)
#' @export dist.xyiStats
#' @author Kristian K Ullrich
#dist.xyiStats<-function(dIUPAC, x.pos, y.pos, i.pos,
#  x.name="x", y.name="y", i.name="i"){
dist.xyiStats<-function(tmpSEQ, x.pos, y.pos, i.pos,
  x.name="x", y.name="y", i.name="i", dist="IUPAC", ncores=1){
    options(scipen=22)
    if(dist=="IUPAC"){
        dIUPAC<-rcpp_distIUPAC(as.character(tmpSEQ), ncores=ncores)
    }
    if(dist!="IUPAC"){
        dIUPAC<-setNames( list(
          dist.dna(as.DNAbin(DNAMultipleAlignment(tmpSEQ)),
            model=dist, as.matrix=TRUE, pairwise.deletion=TRUE),
          pairwiseDeletion(as.character(tmpSEQ))$sitesUsed
          ), c("distIUPAC", "sitesUsed") )
    }
    XNAME<-x.name
    YNAME<-y.name
    INAME<-i.name
    #x
    dMean.x<-NA
    dSd.x<-NA
    dSites.x<-NA
    dMin.x<-NA
    dMax.x<-NA
    #y
    dMean.y<-NA
    dSd.y<-NA
    dSites.y<-NA
    dMin.y<-NA
    dMax.y<-NA
    #i
    dMean.i<-NA
    dSd.i<-NA
    dSites.i<-NA
    dMin.i<-NA
    dMax.i<-NA
    #xy
    dMean.xy<-NA
    dMean_.xy<-NA
    dSd.xy<-NA
    dSites.xy<-NA
    dMin.xy<-NA
    dMax.xy<-NA
    dTotal.xy<-NA
    dSweighted.xy<-NA
    Fst.xy<-NA
    dRelative.xy<-NA
    #xi
    dMean.xi<-NA
    dMean_.xi<-NA
    dSd.xi<-NA
    dSites.xi<-NA
    dMin.xi<-NA
    dMax.xi<-NA
    dTotal.xi<-NA
    dSweighted.xi<-NA
    Fst.xi<-NA
    dRelative.xi<-NA
    #yi
    dMean.yi<-NA
    dMean_.yi<-NA
    dSd.yi<-NA
    dSites.yi<-NA
    dMin.yi<-NA
    dMax.yi<-NA
    dTotal.yi<-NA
    dSweighted.yi<-NA
    Fst.yi<-NA
    dRelative.yi<-NA
    #xyi
    dMean.xyi<-NA
    dSd.xyi<-NA
    dSites.xyi<-NA
    dMin.xyi<-NA
    dMax.xyi<-NA
    dTotal.xyi<-NA
    #xy_xi
    deltaMean.xy_xi<-NA
    deltaSum.xy_xi<-NA
    deltaRatio.xy_xi<-NA
    deltaMean_.xy_xi<-NA
    deltaSum_.xy_xi<-NA
    deltaRatio_.xy_xi<-NA
    deltaMin.xy_xi<-NA
    deltaMinSum.xy_xi<-NA
    deltaMinRatio.xy_xi<-NA
    #yi_xy
    deltaMean.yi_xy<-NA
    deltaSum.yi_xy<-NA
    deltaRatio.yi_xy<-NA
    deltaMean_.yi_xy<-NA
    deltaSum_.yi_xy<-NA
    deltaRatio_.yi_xy<-NA
    deltaMin.yi_xy<-NA
    deltaMinSum.yi_xy<-NA
    deltaMinRatio.yi_xy<-NA
    #yi_xi
    deltaMean.yi_xi<-NA
    deltaSum.yi_xi<-NA
    deltaRatio.yi_xi<-NA
    deltaMean_.yi_xi<-NA
    deltaSum_.yi_xi<-NA
    deltaRatio_.yi_xi<-NA
    deltaMin.yi_xi<-NA
    deltaMinSum.yi_xi<-NA
    deltaMinRatio.yi_xi<-NA
    #
    RND.xyi<-NA
    Gmin.xyi<-NA
    RNDmin.xyi<-NA
    OUT<-list(XNAME, YNAME, INAME,
      dMean.x, dSd.x, dSites.x, dMin.x, dMax.x,
      dMean.y, dSd.y, dSites.y, dMin.y, dMax.y,
      dMean.i, dSd.i, dSites.i, dMin.i, dMax.i,
      dMean.xy, dMean_.xy, dSd.xy, dSites.xy, dMin.xy, dMax.xy,
      dTotal.xy, dSweighted.xy, Fst.xy, dRelative.xy,
      dMean.xi, dMean_.xi, dSd.xi, dSites.xi, dMin.xi, dMax.xi,
      dTotal.xi, dSweighted.xi, Fst.xi, dRelative.xi,
      dMean.yi, dMean_.yi, dSd.yi, dSites.yi, dMin.yi, dMax.yi,
      dTotal.yi, dSweighted.yi, Fst.yi, dRelative.yi,
      dMean.xyi, dSd.xyi, dSites.xyi, dMin.xyi, dMax.xyi, dTotal.xyi,
      deltaMean.xy_xi, deltaSum.xy_xi, deltaRatio.xy_xi,
      deltaMean_.xy_xi, deltaSum_.xy_xi, deltaRatio_.xy_xi,
      deltaMin.xy_xi, deltaMinSum.xy_xi, deltaMinRatio.xy_xi,
      deltaMean.yi_xy, deltaSum.yi_xy, deltaRatio.yi_xy,
      deltaMean_.yi_xy, deltaSum_.yi_xy, deltaRatio_.yi_xy,
      deltaMin.yi_xy, deltaMinSum.yi_xy, deltaMinRatio.yi_xy,
      deltaMean.yi_xi, deltaSum.yi_xi, deltaRatio.yi_xi,
      deltaMean_.yi_xi, deltaSum_.yi_xi, deltaRatio_.yi_xi,
      deltaMin.yi_xi, deltaMinSum.yi_xi, deltaMinRatio.yi_xi,
      RND.xyi, Gmin.xyi, RNDmin.xyi)
    names(OUT)<-c("XNAME", "YNAME", "ONAME",
      "dMean.x", "dSd.x", "dSites.x", "dMin.x", "dMax.x",
      "dMean.y", "dSd.y", "dSites.y", "dMin.y", "dMax.y",
      "dMean.i", "dSd.i", "dSites.i", "dMin.i", "dMax.i",
      "dMean.xy", "dMean_.xy", "dSd.xy", "dSites.xy", "dMin.xy", "dMax.xy",
      "dTotal.xy", "dSweighted.xy", "Fst.xy", "dRelative.xy",
      "dMean.xi", "dMean_.xi", "dSd.xi", "dSites.xi", "dMin.xi", "dMax.xi",
      "dTotal.xi", "dSweighted.xi", "Fst.xi", "dRelative.xi",
      "dMean.yi", "dMean_.yi", "dSd.yi", "dSites.yi", "dMin.yi", "dMax.yi",
      "dTotal.yi", "dSweighted.yi", "Fst.yi", "dRelative.yi",
      "dMean.xyi", "dSd.xyi", "dSites.xyi", "dMin.xyi", "dMax.xyi",
      "dTotal.xyi", "deltaMean.xy_xi", "deltaSum.xy_xi", "deltaRatio.xy_xi",
      "deltaMean_.xy_xi", "deltaSum_.xy_xi", "deltaRatio_.xy_xi",
      "deltaMin.xy_xi", "deltaMinSum.xy_xi", "deltaMinRatio.xy_xi",
      "deltaMean.yi_xy", "deltaSum.yi_xy", "deltaRatio.yi_xy",
      "deltaMean_.yi_xy", "deltaSum_.yi_xy", "deltaRatio_.yi_xy",
      "deltaMin.yi_xy", "deltaMinSum.yi_xy", "deltaMinRatio.yi_xy",
      "deltaMean.yi_xi", "deltaSum.yi_xi", "deltaRatio.yi_xi",
      "deltaMean_.yi_xi", "deltaSum_.yi_xi", "deltaRatio_.yi_xi",
      "deltaMin.yi_xi", "deltaMinSum.yi_xi", "deltaMinRatio.yi_xi",
      "RND.xyi", "Gmin.xyi", "RNDmin.xyi")
    #x
    if(length(x.pos)==1){
        OUT$dMean.x<-mean(dIUPAC$distIUPAC[x.pos, x.pos], na.rm=TRUE)
        OUT$dSd.x<-sd(dIUPAC$distIUPAC[x.pos, x.pos], na.rm=TRUE)
        OUT$dSites.x<-mean(dIUPAC$sitesUsed[x.pos, x.pos], na.rm=TRUE)
        OUT$dMin.x<-min(dIUPAC$distIUPAC[x.pos, x.pos], na.rm=TRUE)
        OUT$dMax.x<-max(dIUPAC$distIUPAC[x.pos, x.pos], na.rm=TRUE)
    }
    else{
        OUT$dMean.x<-mean(as.dist(dIUPAC$distIUPAC[x.pos, x.pos]), na.rm=TRUE)
        OUT$dSd.x<-sd(as.dist(dIUPAC$distIUPAC[x.pos, x.pos]), na.rm=TRUE)
        OUT$dSites.x<-mean(as.dist(dIUPAC$sitesUsed[x.pos, x.pos]), na.rm=TRUE)
        OUT$dMin.x<-min(as.dist(dIUPAC$distIUPAC[x.pos, x.pos]), na.rm=TRUE)
        OUT$dMax.x<-max(as.dist(dIUPAC$distIUPAC[x.pos, x.pos]), na.rm=TRUE)
    }
    #y
    if(length(y.pos)==1){
        OUT$dMean.y<-mean(dIUPAC$distIUPAC[y.pos, y.pos], na.rm=TRUE)
        OUT$dSd.y<-sd(dIUPAC$distIUPAC[y.pos, y.pos], na.rm=TRUE)
        OUT$dSites.y<-mean(dIUPAC$sitesUsed[y.pos, y.pos], na.rm=TRUE)
        OUT$dMin.y<-min(dIUPAC$distIUPAC[y.pos, y.pos], na.rm=TRUE)
        OUT$dMax.y<-max(dIUPAC$distIUPAC[y.pos, y.pos], na.rm=TRUE)
    }
    else{
        OUT$dMean.y<-mean(as.dist(dIUPAC$distIUPAC[y.pos, y.pos]), na.rm=TRUE)
        OUT$dSd.y<-sd(as.dist(dIUPAC$distIUPAC[y.pos, y.pos]), na.rm=TRUE)
        OUT$dSites.y<-mean(as.dist(dIUPAC$sitesUsed[y.pos, y.pos]), na.rm=TRUE)
        OUT$dMin.y<-min(as.dist(dIUPAC$distIUPAC[y.pos, y.pos]), na.rm=TRUE)
        OUT$dMax.y<-max(as.dist(dIUPAC$distIUPAC[y.pos, y.pos]), na.rm=TRUE)
    }
    #i
    if(length(i.pos)==1){
        OUT$dMean.i<-mean(dIUPAC$distIUPAC[i.pos, i.pos], na.rm=TRUE)
        OUT$dSd.i<-sd(dIUPAC$distIUPAC[i.pos, i.pos], na.rm=TRUE)
        OUT$dSites.i<-mean(dIUPAC$sitesUsed[i.pos, i.pos], na.rm=TRUE)
        OUT$dMin.i<-min(dIUPAC$distIUPAC[i.pos, i.pos], na.rm=TRUE)
        OUT$dMax.i<-max(dIUPAC$distIUPAC[i.pos, i.pos], na.rm=TRUE)
    }
    else{
        OUT$dMean.i<-mean(as.dist(dIUPAC$distIUPAC[i.pos, i.pos]), na.rm=TRUE)
        OUT$dSd.i<-sd(as.dist(dIUPAC$distIUPAC[i.pos, i.pos]), na.rm=TRUE)
        OUT$dSites.i<-mean(as.dist(dIUPAC$sitesUsed[i.pos, i.pos]), na.rm=TRUE)
        OUT$dMin.i<-min(as.dist(dIUPAC$distIUPAC[i.pos, i.pos]), na.rm=TRUE)
        OUT$dMax.i<-max(as.dist(dIUPAC$distIUPAC[i.pos, i.pos]), na.rm=TRUE)
    }
    #xy
    OUT$dMean.xy<-mean(dIUPAC$distIUPAC[x.pos, y.pos], na.rm=TRUE)
    OUT$dMean_.xy<-OUT$dMean.xy - (OUT$dMean.x/2) - (OUT$dMean.y/2) #F2(x,y) equation (17) Peter, Genetics, Vol. 202, 1485-1501, April 2016
    OUT$dSd.xy<-sd(dIUPAC$distIUPAC[x.pos, y.pos], na.rm=TRUE)
    OUT$dSites.xy<-mean(dIUPAC$sitesUsed[x.pos, y.pos], na.rm=TRUE)
    OUT$dMin.xy<-min(dIUPAC$distIUPAC[x.pos, y.pos], na.rm=TRUE)
    OUT$dMax.xy<-max(dIUPAC$distIUPAC[x.pos, y.pos], na.rm=TRUE)
    OUT$dTotal.xy<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, y.pos),
      c(x.pos, y.pos)]), na.rm=TRUE)
    OUT$dSweighted.xy<-((length(x.pos)/(length(c(x.pos, y.pos)))) * 
      OUT$dMean.x) + ((length(y.pos)/(length(c(x.pos, y.pos)))) * OUT$dMean.y)
    OUT$Fst.xy<-(OUT$dTotal.xy - OUT$dSweighted.xy) / OUT$dTotal.xy
    OUT$dRelative.xy<-OUT$dMean.xy - OUT$dSweighted.xy
    #xi
    OUT$dMean.xi<-mean(dIUPAC$distIUPAC[x.pos, i.pos], na.rm=TRUE)
    OUT$dMean_.xi<-OUT$dMean.xi - (OUT$dMean.x/2) - (OUT$dMean.i/2) #F2(x,i) equation (17) Peter, Genetics, Vol. 202, 1485-1501, April 2016
    OUT$dSd.xi<-sd(dIUPAC$distIUPAC[x.pos, i.pos], na.rm=TRUE)
    OUT$dSites.xi<-mean(dIUPAC$sitesUsed[x.pos, i.pos], na.rm=TRUE)
    OUT$dMin.xi<-min(dIUPAC$distIUPAC[x.pos, i.pos], na.rm=TRUE)
    OUT$dMax.xi<-max(dIUPAC$distIUPAC[x.pos, i.pos], na.rm=TRUE)
    OUT$dTotal.xi<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, i.pos),
      c(x.pos, i.pos)]), na.rm=TRUE)
    OUT$dSweighted.xi<-((length(x.pos)/(length(c(x.pos, i.pos)))) *
      OUT$dMean.x) + ((length(i.pos)/(length(c(x.pos, i.pos)))) * OUT$dMean.i)
    OUT$Fst.xi<-(OUT$dTotal.xi - OUT$dSweighted.xi) / OUT$dTotal.xi
    OUT$dRelative.xi<-OUT$dMean.xi - OUT$dSweighted.xi
    #yi
    OUT$dMean.yi<-mean(dIUPAC$distIUPAC[y.pos, i.pos], na.rm=TRUE)
    OUT$dMean_.yi<-OUT$dMean.yi - (OUT$dMean.y/2) - (OUT$dMean.i/2) #F2(y,i) equation (17) Peter, Genetics, Vol. 202, 1485-1501, April 2016
    OUT$dSd.yi<-sd(dIUPAC$distIUPAC[y.pos, i.pos], na.rm=TRUE)
    OUT$dSites.yi<-mean(dIUPAC$sitesUsed[y.pos, i.pos], na.rm=TRUE)
    OUT$dMin.yi<-min(dIUPAC$distIUPAC[y.pos, i.pos], na.rm=TRUE)
    OUT$dMax.yi<-max(dIUPAC$distIUPAC[y.pos, i.pos], na.rm=TRUE)
    OUT$dTotal.yi<-mean(as.dist(dIUPAC$distIUPAC[c(y.pos, i.pos),
      c(y.pos, i.pos)]), na.rm=TRUE)
    OUT$dSweighted.yi<-((length(y.pos)/(length(c(y.pos, i.pos)))) *
      OUT$dMean.y) + ((length(i.pos)/(length(c(y.pos, i.pos)))) * OUT$dMean.i)
    OUT$Fst.yi<-(OUT$dTotal.yi - OUT$dSweighted.yi) / OUT$dTotal.yi
    OUT$dRelative.yi<-OUT$dMean.yi - OUT$dSweighted.yi
    #xyi
    OUT$dMean.xyi<-mean(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos)],
      dIUPAC$distIUPAC[y.pos, i.pos]), na.rm=TRUE)
    OUT$dSd.xyi<-sd(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos)],
      dIUPAC$distIUPAC[y.pos, i.pos]), na.rm=TRUE)
    OUT$dSites.xyi<-mean(c(dIUPAC$sitesUsed[x.pos, c(y.pos, i.pos)],
      dIUPAC$sitesUsed[y.pos, i.pos]), na.rm=TRUE)
    OUT$dMin.xyi<-min(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos)],
      dIUPAC$distIUPAC[y.pos, i.pos]), na.rm=TRUE)
    OUT$dMax.xyi<-max(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos)],
      dIUPAC$distIUPAC[y.pos, i.pos]), na.rm=TRUE)
    OUT$dTotal.xyi<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, y.pos, i.pos), 
      c(x.pos, y.pos, i.pos)]), na.rm=TRUE)
    #xy_xi
    OUT$deltaMean.xy_xi<-OUT$dMean.xy - OUT$dMean.xi
    OUT$deltaSum.xy_xi<-OUT$dMean.xy + OUT$dMean.xi
    OUT$deltaRatio.xy_xi<-OUT$deltaMean.xy_xi / OUT$deltaSum.xy_xi
    OUT$deltaMean_.xy_xi<-OUT$dMean_.xy - OUT$dMean_.xi
    OUT$deltaSum_.xy_xi<-OUT$dMean_.xy + OUT$dMean_.xi
    OUT$deltaRatio_.xy_xi<-OUT$deltaMean_.xy_xi / OUT$deltaSum_.xy_xi
    OUT$deltaMin.xy_xi<-OUT$dMin.xy - OUT$dMean.xi
    OUT$deltaMinSum.xy_xi<-OUT$dMin.xy + OUT$dMean.xi
    OUT$deltaMinRatio.xy_xi<-OUT$deltaMin.xy_xi / OUT$deltaMinSum.xy_xi
    #yi_xy
    OUT$deltaMean.yi_xy<-OUT$dMean.yi - OUT$dMean.xy
    OUT$deltaSum.yi_xy<-OUT$dMean.yi + OUT$dMean.xy
    OUT$deltaRatio.yi_xy<-OUT$deltaMean.yi_xy / OUT$deltaSum.yi_xy
    OUT$deltaMean_.yi_xy<-OUT$dMean_.yi - OUT$dMean_.xy
    OUT$deltaSum_.yi_xy<-OUT$dMean_.yi + OUT$dMean_.xy
    OUT$deltaRatio_.yi_xy<-OUT$deltaMean_.yi_xy / OUT$deltaSum_.yi_xy
    OUT$deltaMin.yi_xy<-OUT$dMin.yi - OUT$dMean.xy
    OUT$deltaMinSum.yi_xy<-OUT$dMin.yi + OUT$dMean.xy
    OUT$deltaMinRatio.yi_xy<-OUT$deltaMin.yi_xy / OUT$deltaMinSum.yi_xy
    #yi_xi
    OUT$deltaMean.yi_xi<-OUT$dMean.yi - OUT$dMean.xi
    OUT$deltaSum.yi_xi<-OUT$dMean.yi + OUT$dMean.xi
    OUT$deltaRatio.yi_xi<-OUT$deltaMean.yi_xi / OUT$deltaSum.yi_xi
    OUT$deltaMean_.yi_xi<-OUT$dMean_.yi - OUT$dMean_.xi
    OUT$deltaSum_.yi_xi<-OUT$dMean_.yi + OUT$dMean_.xi
    OUT$deltaRatio_.yi_xi<-OUT$deltaMean_.yi_xi / OUT$deltaSum_.yi_xi
    OUT$deltaMin.yi_xi<-OUT$dMin.yi - OUT$dMean.xi
    OUT$deltaMinSum.yi_xi<-OUT$dMin.yi + OUT$dMean.xi
    OUT$deltaMinRatio.yi_xi<-OUT$deltaMin.yi_xi / OUT$deltaMinSum.yi_xi
    #
    OUT$RND.xyi<-OUT$dMean.xy/((OUT$dMean.xi + OUT$dMean.yi)/2)
    OUT$Gmin.xyi<-OUT$dMin.xy/OUT$dMean.xy
    OUT$RNDmin.xyi<-OUT$dMin.xy/((OUT$dMean.xi + OUT$dMean.yi)/2)
    return(OUT)
}
