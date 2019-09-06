#' @title dist.xyoStats
#' @name dist.xyoStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver; y: donor)
#' with an ingroup population (i: ingroup).
#' In a four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P2 y:P3 i:O].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow
#' from P2>>P3, the populations should be defined
#' for RND, Gmin and RNDmin statistics as follows [x:P3 y:P2 i:O].
#' @importFrom stats as.dist sd setNames
#' @param tmpSEQ \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [mandatory]
#' @param y.pos population Y positions [mandatory]
#' @param o.pos population I positions [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param o.name population I name [default: "o"]
#' @param dist distance to use [default: IUPAC] or choose one model as in
#' \link[ape]{dist.dna} [default: "IUPAC"]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC} [default: 1]
#' @seealso \code{\link[distIUPAC]{xyoStats}},
#' \code{\link[distIUPAC]{distIUPAC}}, \link[distIUPAC]{rcpp_distIUPAC}
#' @examples
#' ##Use the 'xyoStats' function to handle population assignment automatically
#' 
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' 
#' ##Here, one needs to consider the changed x, y and o positions due
#' ##to sub-sampling.
#' CAS.pos_<-seq_along(CAS.pos)
#' AFG.pos_<-seq(from=length(CAS.pos_)+1, to=length(c(CAS.pos, AFG.pos)))
#' SPRE.pos_<-seq(from=length(c(CAS.pos, AFG.pos))+1,
#'   to=length(c(CAS.pos, AFG.pos, SPRE.pos)))
#' 
#' ##pairwise deletion
#' CAS_AFG_SPRE.pairwiseDeletion.dist<-distIUPAC(
#'   as.character(MySequences[c(CAS.pos,AFG.pos,SPRE.pos)]))
#' CAS_AFG_SPRE.pairwiseDeletion.xyoStats<-dist.xyoStats(
#'   CAS_AFG_SPRE.pairwiseDeletion.dist, 
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, o.pos=SPRE.pos_,
#'   x.name="CAS", y.name="AFG", o.name="SPRE")
#' 
#' ##global deletion
#' CAS_AFG_SPRE.globalDeletion.dist<-distIUPAC(
#'   as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos,SPRE.pos)])))
#' CAS_AFG_SPRE.globalDeletion.xyoStats<-dist.xyoStats(
#'   CAS_AFG_SPRE.globalDeletion.dist,
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, o.pos=SPRE.pos_,
#'   x.name="CAS", y.name="AFG", o.name="SPRE")
#' 
#' ##compare results
#' rbind(CAS_AFG_SPRE.pairwiseDeletion.xyoStats,
#'   CAS_AFG_SPRE.globalDeletion.xyoStats)
#' @export dist.xyoStats
#' @author Kristian K Ullrich
#dist.xyoStats<-function(dIUPAC, x.pos, y.pos, o.pos,
#  x.name="x", y.name="y", o.name="o"){
dist.xyoStats<-function(tmpSEQ, x.pos, y.pos, o.pos,
  x.name="x", y.name="y", o.name="o", dist="IUPAC", ncores=1){
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
    ONAME<-o.name
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
    #o
    dMean.o<-NA
    dSd.o<-NA
    dSites.o<-NA
    dMin.o<-NA
    dMax.o<-NA
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
    #xo
    dMean.xo<-NA
    dMean_.xo<-NA
    dSd.xo<-NA
    dSites.xo<-NA
    dMin.xo<-NA
    dMax.xo<-NA
    dTotal.xo<-NA
    dSweighted.xo<-NA
    Fst.xo<-NA
    dRelative.xo<-NA
    #yo
    dMean.yo<-NA
    dMean_.yo<-NA
    dSd.yo<-NA
    dSites.yo<-NA
    dMin.yo<-NA
    dMax.yo<-NA
    dTotal.yo<-NA
    dSweighted.yo<-NA
    Fst.yo<-NA
    dRelative.yo<-NA
    #xyo
    dMean.xyo<-NA
    dSd.xyo<-NA
    dSites.xyo<-NA
    dMin.xyo<-NA
    dMax.xyo<-NA
    dTotal.xyo<-NA
    #xy_xo
    deltaMean.xy_xo<-NA
    deltaSum.xy_xo<-NA
    deltaRatio.xy_xo<-NA
    deltaMean_.xy_xo<-NA
    deltaSum_.xy_xo<-NA
    deltaRatio_.xy_xo<-NA
    deltaMin.xy_xo<-NA
    deltaMinSum.xy_xo<-NA
    deltaMinRatio.xy_xo<-NA
    #yo_xy
    deltaMean.yo_xy<-NA
    deltaSum.yo_xy<-NA
    deltaRatio.yo_xy<-NA
    deltaMean_.yo_xy<-NA
    deltaSum_.yo_xy<-NA
    deltaRatio_.yo_xy<-NA
    deltaMin.yo_xy<-NA
    deltaMinSum.yo_xy<-NA
    deltaMinRatio.yo_xy<-NA
    #yo_xo
    deltaMean.yo_xo<-NA
    deltaSum.yo_xo<-NA
    deltaRatio.yo_xo<-NA
    deltaMean_.yo_xo<-NA
    deltaSum_.yo_xo<-NA
    deltaRatio_.yo_xo<-NA
    deltaMin.yo_xo<-NA
    deltaMinSum.yo_xo<-NA
    deltaMinRatio.yo_xo<-NA
    #
    RND.xyo<-NA
    Gmin.xyo<-NA
    RNDmin.xyo<-NA
    OUT<-list(XNAME, YNAME, ONAME,
      dMean.x, dSd.x, dSites.x, dMin.x, dMax.x,
      dMean.y, dSd.y, dSites.y, dMin.y, dMax.y,
      dMean.o, dSd.o, dSites.o, dMin.o, dMax.o,
      dMean.xy, dMean_.xy, dSd.xy, dSites.xy, dMin.xy, dMax.xy,
      dTotal.xy, dSweighted.xy, Fst.xy, dRelative.xy,
      dMean.xo, dMean_.xo, dSd.xo, dSites.xo, dMin.xo, dMax.xo,
      dTotal.xo, dSweighted.xo, Fst.xo, dRelative.xo,
      dMean.yo, dMean_.yo, dSd.yo, dSites.yo, dMin.yo, dMax.yo,
      dTotal.yo, dSweighted.yo, Fst.yo, dRelative.yo,
      dMean.xyo, dSd.xyo, dSites.xyo, dMin.xyo, dMax.xyo, dTotal.xyo,
      deltaMean.xy_xo, deltaSum.xy_xo, deltaRatio.xy_xo,
      deltaMean_.xy_xo, deltaSum_.xy_xo, deltaRatio_.xy_xo,
      deltaMin.xy_xo, deltaMinSum.xy_xo, deltaMinRatio.xy_xo,
      deltaMean.yo_xy, deltaSum.yo_xy, deltaRatio.yo_xy,
      deltaMean_.yo_xy, deltaSum_.yo_xy, deltaRatio_.yo_xy,
      deltaMin.yo_xy, deltaMinSum.yo_xy, deltaMinRatio.yo_xy,
      deltaMean.yo_xo, deltaSum.yo_xo, deltaRatio.yo_xo,
      deltaMean_.yo_xo, deltaSum_.yo_xo, deltaRatio_.yo_xo,
      deltaMin.yo_xo, deltaMinSum.yo_xo, deltaMinRatio.yo_xo,
      RND.xyo, Gmin.xyo, RNDmin.xyo)
    names(OUT)<-c("XNAME", "YNAME", "ONAME",
      "dMean.x", "dSd.x", "dSites.x", "dMin.x", "dMax.x",
      "dMean.y", "dSd.y", "dSites.y", "dMin.y", "dMax.y",
      "dMean.o", "dSd.o", "dSites.o", "dMin.o", "dMax.o",
      "dMean.xy", "dMean_.xy", "dSd.xy", "dSites.xy", "dMin.xy", "dMax.xy",
      "dTotal.xy", "dSweighted.xy", "Fst.xy", "dRelative.xy",
      "dMean.xo", "dMean_.xo", "dSd.xo", "dSites.xo", "dMin.xo", "dMax.xo",
      "dTotal.xo", "dSweighted.xo", "Fst.xo", "dRelative.xo",
      "dMean.yo", "dMean_.yo", "dSd.yo", "dSites.yo", "dMin.yo", "dMax.yo",
      "dTotal.yo", "dSweighted.yo", "Fst.yo", "dRelative.yo",
      "dMean.xyo", "dSd.xyo", "dSites.xyo", "dMin.xyo", "dMax.xyo",
      "dTotal.xyo", "deltaMean.xy_xo", "deltaSum.xy_xo", "deltaRatio.xy_xo",
      "deltaMean_.xy_xo", "deltaSum_.xy_xo", "deltaRatio_.xy_xo",
      "deltaMin.xy_xo", "deltaMinSum.xy_xo", "deltaMinRatio.xy_xo",
      "deltaMean.yo_xy", "deltaSum.yo_xy", "deltaRatio.yo_xy",
      "deltaMean_.yo_xy", "deltaSum_.yo_xy", "deltaRatio_.yo_xy",
      "deltaMin.yo_xy", "deltaMinSum.yo_xy", "deltaMinRatio.yo_xy",
      "deltaMean.yo_xo", "deltaSum.yo_xo", "deltaRatio.yo_xo",
      "deltaMean_.yo_xo", "deltaSum_.yo_xo", "deltaRatio_.yo_xo",
      "deltaMin.yo_xo", "deltaMinSum.yo_xo", "deltaMinRatio.yo_xo",
      "RND.xyo","Gmin.xyo","RNDmin.xyo")
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
    #o
    if(length(o.pos)==1){
        OUT$dMean.o<-mean(dIUPAC$distIUPAC[o.pos, o.pos], na.rm=TRUE)
        OUT$dSd.o<-sd(dIUPAC$distIUPAC[o.pos, o.pos], na.rm=TRUE)
        OUT$dSites.o<-mean(dIUPAC$sitesUsed[o.pos, o.pos], na.rm=TRUE)
        OUT$dMin.o<-min(dIUPAC$distIUPAC[o.pos, o.pos], na.rm=TRUE)
        OUT$dMax.o<-max(dIUPAC$distIUPAC[o.pos, o.pos], na.rm=TRUE)
    }
    else{
        OUT$dMean.o<-mean(as.dist(dIUPAC$distIUPAC[o.pos, o.pos]), na.rm=TRUE)
        OUT$dSd.o<-sd(as.dist(dIUPAC$distIUPAC[o.pos, o.pos]), na.rm=TRUE)
        OUT$dSites.o<-mean(as.dist(dIUPAC$sitesUsed[o.pos, o.pos]), na.rm=TRUE)
        OUT$dMin.o<-min(as.dist(dIUPAC$distIUPAC[o.pos, o.pos]), na.rm=TRUE)
        OUT$dMax.o<-max(as.dist(dIUPAC$distIUPAC[o.pos, o.pos]), na.rm=TRUE)
    }
    #xy
    OUT$dMean.xy<-mean(dIUPAC$distIUPAC[x.pos, y.pos], na.rm=TRUE)
    OUT$dMean_.xy<-OUT$dMean.xy - (OUT$dMean.x/2) - (OUT$dMean.y/2)
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
    #xo
    OUT$dMean.xo<-mean(dIUPAC$distIUPAC[x.pos, o.pos], na.rm=TRUE)
    OUT$dMean_.xo<-OUT$dMean.xo - (OUT$dMean.x/2) - (OUT$dMean.o/2)
    OUT$dSd.xo<-sd(dIUPAC$distIUPAC[x.pos, o.pos], na.rm=TRUE)
    OUT$dSites.xo<-mean(dIUPAC$sitesUsed[x.pos, o.pos], na.rm=TRUE)
    OUT$dMin.xo<-min(dIUPAC$distIUPAC[x.pos, o.pos], na.rm=TRUE)
    OUT$dMax.xo<-max(dIUPAC$distIUPAC[x.pos, o.pos], na.rm=TRUE)
    OUT$dTotal.xo<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, o.pos),
      c(x.pos, o.pos)]), na.rm=TRUE)
    OUT$dSweighted.xo<-((length(x.pos)/(length(c(x.pos, o.pos)))) *
      OUT$dMean.x) + ((length(o.pos)/(length(c(x.pos, o.pos)))) * OUT$dMean.o)
    OUT$Fst.xo<-(OUT$dTotal.xo - OUT$dSweighted.xo) / OUT$dTotal.xo
    OUT$dRelative.xo<-OUT$dMean.xo - OUT$dSweighted.xo
    #yo
    OUT$dMean.yo<-mean(dIUPAC$distIUPAC[y.pos, o.pos], na.rm=TRUE)
    OUT$dMean_.yo<-OUT$dMean.yo - (OUT$dMean.y/2) - (OUT$dMean.o/2)
    OUT$dSd.yo<-sd(dIUPAC$distIUPAC[y.pos, o.pos], na.rm=TRUE)
    OUT$dSites.yo<-mean(dIUPAC$sitesUsed[y.pos, o.pos], na.rm=TRUE)
    OUT$dMin.yo<-min(dIUPAC$distIUPAC[y.pos, o.pos], na.rm=TRUE)
    OUT$dMax.yo<-max(dIUPAC$distIUPAC[y.pos, o.pos], na.rm=TRUE)
    OUT$dTotal.yo<-mean(as.dist(dIUPAC$distIUPAC[c(y.pos, o.pos),
      c(y.pos, o.pos)]), na.rm=TRUE)
    OUT$dSweighted.yo<-((length(y.pos)/(length(c(y.pos, o.pos)))) *
      OUT$dMean.y) + ((length(o.pos)/(length(c(y.pos, o.pos)))) * OUT$dMean.o)
    OUT$Fst.yo<-(OUT$dTotal.yo - OUT$dSweighted.yo) / OUT$dTotal.yo
    OUT$dRelative.yo<-OUT$dMean.yo - OUT$dSweighted.yo
    #xyo
    OUT$dMean.xyo<-mean(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dSd.xyo<-sd(c(dIUPAC$distIUPAC[x.pos,c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dSites.xyo<-mean(c(dIUPAC$sitesUsed[x.pos, c(y.pos, o.pos)],
      dIUPAC$sitesUsed[y.pos, o.pos]), na.rm=TRUE)
    OUT$dMin.xyo<-min(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dMax.xyo<-max(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dTotal.xyo<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, y.pos,o.pos),
      c(x.pos, y.pos,o.pos)]), na.rm=TRUE)
    #xy_xo
    OUT$deltaMean.xy_xo<-OUT$dMean.xy - OUT$dMean.xo
    OUT$deltaSum.xy_xo<-OUT$dMean.xy + OUT$dMean.xo
    OUT$deltaRatio.xy_xo<-OUT$deltaMean.xy_xo / OUT$deltaSum.xy_xo
    OUT$deltaMean_.xy_xo<-OUT$dMean_.xy - OUT$dMean_.xo
    OUT$deltaSum_.xy_xo<-OUT$dMean_.xy + OUT$dMean_.xo
    OUT$deltaRatio_.xy_xo<-OUT$deltaMean_.xy_xo / OUT$deltaSum_.xy_xo
    OUT$deltaMin.xy_xo<-OUT$dMin.xy - OUT$dMean.xo
    OUT$deltaMinSum.xy_xo<-OUT$dMin.xy + OUT$dMean.xo
    OUT$deltaMinRatio.xy_xo<-OUT$deltaMin.xy_xo / OUT$deltaMinSum.xy_xo
    #yo_xy
    OUT$deltaMean.yo_xy<-OUT$dMean.yo - OUT$dMean.xy
    OUT$deltaSum.yo_xy<-OUT$dMean.yo + OUT$dMean.xy
    OUT$deltaRatio.yo_xy<-OUT$deltaMean.yo_xy / OUT$deltaSum.yo_xy
    OUT$deltaMean_.yo_xy<-OUT$dMean_.yo - OUT$dMean_.xy
    OUT$deltaSum_.yo_xy<-OUT$dMean_.yo + OUT$dMean_.xy
    OUT$deltaRatio_.yo_xy<-OUT$deltaMean_.yo_xy / OUT$deltaSum_.yo_xy
    OUT$deltaMin.yo_xy<-OUT$dMin.yo - OUT$dMean.xy
    OUT$deltaMinSum.yo_xy<-OUT$dMin.yo + OUT$dMean.xy
    OUT$deltaMinRatio.yo_xy<-OUT$deltaMin.yo_xy / OUT$deltaMinSum.yo_xy
    #yo_xo
    OUT$deltaMean.yo_xo<-OUT$dMean.yo - OUT$dMean.xo
    OUT$deltaSum.yo_xo<-OUT$dMean.yo + OUT$dMean.xo
    OUT$deltaRatio.yo_xo<-OUT$deltaMean.yo_xo / OUT$deltaSum.yo_xo
    OUT$deltaMean_.yo_xo<-OUT$dMean_.yo - OUT$dMean_.xo
    OUT$deltaSum_.yo_xo<-OUT$dMean_.yo + OUT$dMean_.xo
    OUT$deltaRatio_.yo_xo<-OUT$deltaMean_.yo_xo / OUT$deltaSum_.yo_xo
    OUT$deltaMin.yo_xo<-OUT$dMin.yo - OUT$dMean.xo
    OUT$deltaMinSum.yo_xo<-OUT$dMin.yo + OUT$dMean.xo
    OUT$deltaMinRatio.yo_xo<-OUT$deltaMin.yo_xo / OUT$deltaMinSum.yo_xo
    #
    OUT$RND.xyo<-OUT$dMean.xy/((OUT$dMean.xo + OUT$dMean.yo)/2)
    OUT$Gmin.xyo<-OUT$dMin.xy/OUT$dMean.xy
    OUT$RNDmin.xyo<-OUT$dMin.xy/((OUT$dMean.xo + OUT$dMean.yo)/2)
    return(OUT)
}
