#' @title abbababa.xyioStats
#' @name abbababa.xyioStats
#' @description This function calculates ABBA-BABA statistics
#' comparing two populations (x: receiver; y: donor) with an ingroup
#' population (i: ingroup) and an outgroup population (o: outgroup).
#' In a four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined as follows [x:P2 y:P3 i:P1 o:P4].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined as follows [x:P3 y:P2 i:P1 o:P4].
#' @importFrom stats as.dist sd setNames
#' @param tmpSEQ \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [mandatory]
#' @param y.pos population Y positions [mandatory]
#' @param i.pos population I positions [mandatory]
#' @param o.pos population O positions [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param i.name population I name [default: "i"]
#' @param o.name population I name [default: "o"]
#' @param dist distance to use [default: IUPAC] or choose one model as in
#' \link[ape]{dist.dna} [default: "IUPAC"]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC} [default: 1]
#' @seealso \code{\link[distIUPAC]{xyioStats}},
#' \code{\link[distIUPAC]{distIUPAC}}, \code{\link[distIUPAC]{rcpp_distIUPAC}}
#' @examples
#' ##Use the 'xyioStats' function to handle population assignment automatically
#' 
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' IRA.pos<-71:78
#' APO.pos<-1
#' 
#' ##Here, one needs to consider the changed x, y and i positions due
#' ##to sub-sampling.
#' CAS.pos_<-seq_along(CAS.pos)
#' AFG.pos_<-seq(from=length(CAS.pos_)+1, to=length(c(CAS.pos, AFG.pos)))
#' IRA.pos_<-seq(from=length(c(CAS.pos, AFG.pos))+1,
#' to=length(c(CAS.pos, AFG.pos, IRA.pos)))
#' APO.pos_<-seq(from=length(c(CAS.pos, AFG.pos, IRA.pos))+1,
#' to=length(c(CAS.pos, AFG.pos, IRA.pos, APO.pos)))
#' 
#' ##pairwise deletion
#' CAS_AFG_IRA_APO.pairwiseDeletion.dist<-distIUPAC(
#'   as.character(MySequences[c(CAS.pos,AFG.pos,IRA.pos,APO.pos)]))
#' CAS_AFG_IRA_APO.pairwiseDeletion.xyioStats<-dist.xyioStats(
#'   CAS_AFG_IRA_APO.pairwiseDeletion.dist, 
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_, o.pos=APO.pos_,
#'   x.name="CAS", y.name="AFG", i.name="IRA", o.name="APO")
#' 
#' ##global deletion
#' CAS_AFG_IRA_APO.globalDeletion.dist<-distIUPAC(
#'   as.character(globalDeletion(MySequences[c(CAS.pos,AFG.pos,
#'   IRA.pos,APO.pos)])))
#' CAS_AFG_IRA_APO.globalDeletion.xyioStats<-dist.xyioStats(
#'   CAS_AFG_IRA_APO.globalDeletion.dist,
#'   x.pos=CAS.pos_, y.pos=AFG.pos_, i.pos=IRA.pos_, o.pos=APO.pos_,
#'   x.name="CAS", y.name="AFG", i.name="IRA", o.name="APO")
#' 
#' ##compare results
#' rbind(CAS_AFG_IRA_APO.pairwiseDeletion.xyiStats,
#'   CAS_AFG_IRA_APO.globalDeletion.xyiStats)
#' @export abbababa.xyioStats
#' @author Kristian K Ullrich
abbababa.xyioStats<-function(tmpSEQ, x.pos, y.pos, i.pos, o.pos,
  x.freq=1.0, y.freq=1.0, i.freq=1.0, o.freq=1.0,
  x.name="x", y.name="y", i.name="i", o.name="o", dist="IUPAC", ncores=1){
    options(scipen=22)
    
    XNAME<-x.name
    YNAME<-y.name
    INAME<-i.name
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
    #i
    dMean.i<-NA
    dSd.i<-NA
    dSites.i<-NA
    dMin.i<-NA
    dMax.i<-NA
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
    #io
    dMean.io<-NA
    dMean_.io<-NA
    dSd.io<-NA
    dSites.io<-NA
    dMin.io<-NA
    dMax.io<-NA
    dTotal.io<-NA
    dSweighted.io<-NA
    Fst.io<-NA
    dRelative.io<-NA
    #xyi
    dMean.xyi<-NA
    dSd.xyi<-NA
    dSites.xyi<-NA
    dMin.xyi<-NA
    dMax.xyi<-NA
    dTotal.xyi<-NA
    #xyo
    dMean.xyo<-NA
    dSd.xyo<-NA
    dSites.xyo<-NA
    dMin.xyo<-NA
    dMax.xyo<-NA
    dTotal.xyo<-NA
    #yio
    dMean.yio<-NA
    dSd.yio<-NA
    dSites.yio<-NA
    dMin.yio<-NA
    dMax.yio<-NA
    dTotal.yio<-NA
    #xyio
    dMean.xyio<-NA
    dSd.xyio<-NA
    dSites.xyio<-NA
    dMin.xyio<-NA
    dMax.xyio<-NA
    dTotal.xyio<-NA
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
    #
    RND.xyo<-NA
    Gmin.xyo<-NA
    RNDmin.xyo<-NA
    OUT<-list(XNAME, YNAME, INAME, ONAME,
      dMean.x, dSd.x, dSites.x, dMin.x, dMax.x,
      dMean.y, dSd.y, dSites.y, dMin.y, dMax.y,
      dMean.i, dSd.i, dSites.i, dMin.i, dMax.i,
      dMean.o, dSd.o, dSites.o, dMin.o, dMax.o,
      #xy
      dMean.xy, dMean_.xy, dSd.xy, dSites.xy, dMin.xy, dMax.xy,
      dTotal.xy, dSweighted.xy, Fst.xy, dRelative.xy,
      #xi
      dMean.xi, dMean_.xi, dSd.xi, dSites.xi, dMin.xi, dMax.xi,
      dTotal.xi, dSweighted.xi, Fst.xi, dRelative.xi,
      #xo
      dMean.xo, dMean_.xo, dSd.xo, dSites.xo, dMin.xo, dMax.xo,
      dTotal.xo, dSweighted.xo, Fst.xo, dRelative.xo,
      #yi
      dMean.yi, dMean_.yi, dSd.yi, dSites.yi, dMin.yi, dMax.yi,
      dTotal.yi, dSweighted.yi, Fst.yi, dRelative.yi,
      #yo
      dMean.yo, dMean_.yo, dSd.yo, dSites.yo, dMin.yo, dMax.yo,
      dTotal.yo, dSweighted.yo, Fst.yo, dRelative.yo,
      #io
      dMean.io, dMean_.io, dSd.io, dSites.io, dMin.io, dMax.io,
      dTotal.io, dSweighted.io, Fst.io, dRelative.io,
      #xyi
      dMean.xyi, dSd.xyi, dSites.xyi, dMin.xyi, dMax.xyi, dTotal.xyi,
      #xyo
      dMean.xyo, dSd.xyo, dSites.xyo, dMin.xyo, dMax.xyo, dTotal.xyo,
      #yio
      dMean.yio, dSd.yio, dSites.yio, dMin.yio, dMax.yio, dTotal.yio,
      #xyio
      dMean.xyio, dSd.xyio, dSites.xyio, dMin.xyio, dMax.xyio, dTotal.xyio,
      #xy_xi
      deltaMean.xy_xi, deltaSum.xy_xi, deltaRatio.xy_xi,
      deltaMean_.xy_xi, deltaSum_.xy_xi, deltaRatio_.xy_xi,
      deltaMin.xy_xi, deltaMinSum.xy_xi, deltaMinRatio.xy_xi,
      #yi_xy
      deltaMean.yi_xy, deltaSum.yi_xy, deltaRatio.yi_xy,
      deltaMean_.yi_xy, deltaSum_.yi_xy, deltaRatio_.yi_xy,
      deltaMin.yi_xy, deltaMinSum.yi_xy, deltaMinRatio.yi_xy,
      #yi_xi
      deltaMean.yi_xi, deltaSum.yi_xi, deltaRatio.yi_xi,
      deltaMean_.yi_xi, deltaSum_.yi_xi, deltaRatio_.yi_xi,
      deltaMin.yi_xi, deltaMinSum.yi_xi, deltaMinRatio.yi_xi,
      #
      RND.xyi, Gmin.xyi, RNDmin.xyi,
      #
      RND.xyo, Gmin.xyo, RNDmin.xyo)
    names(OUT)<-c("XNAME", "YNAME", "INAME", "ONAME",
      "dMean.x", "dSd.x", "dSites.x", "dMin.x", "dMax.x",
      "dMean.y", "dSd.y", "dSites.y", "dMin.y", "dMax.y",
      "dMean.i", "dSd.i", "dSites.i", "dMin.i", "dMax.i",
      "dMean.o", "dSd.o", "dSites.o", "dMin.o", "dMax.o",
      #xy
      "dMean.xy", "dMean_.xy", "dSd.xy", "dSites.xy", "dMin.xy", "dMax.xy",
      "dTotal.xy", "dSweighted.xy", "Fst.xy", "dRelative.xy",
      #xi
      "dMean.xi", "dMean_.xi", "dSd.xi", "dSites.xi", "dMin.xi", "dMax.xi",
      "dTotal.xi", "dSweighted.xi", "Fst.xi", "dRelative.xi",
      #xo
      "dMean.xo", "dMean_.xo", "dSd.xo", "dSites.xo", "dMin.xo", "dMax.xo",
      "dTotal.xo", "dSweighted.xo", "Fst.xo", "dRelative.xo",
      #yi
      "dMean.yi", "dMean_.yi", "dSd.yi", "dSites.yi", "dMin.yi", "dMax.yi",
      "dTotal.yi", "dSweighted.yi", "Fst.yi", "dRelative.yi",
      #yo
      "dMean.yo", "dMean_.yo", "dSd.yo", "dSites.yo", "dMin.yo", "dMax.yo",
      "dTotal.yo", "dSweighted.yo", "Fst.yo", "dRelative.yo",
      #io
      "dMean.io", "dMean_.io", "dSd.io", "dSites.io", "dMin.io", "dMax.io",
      "dTotal.io", "dSweighted.io", "Fst.io", "dRelative.io",
      #xyi
      "dMean.xyi", "dSd.xyi", "dSites.xyi", "dMin.xyi", "dMax.xyi",
      "dTotal.xyi",
      #xyo
      "dMean.xyo", "dSd.xyo", "dSites.xyo", "dMin.xyo", "dMax.xyo",
      "dTotal.xyo",
      #yio
      "dMean.yio", "dSd.yio", "dSites.yio", "dMin.yio", "dMax.yio",
      "dTotal.yio",
      #xyio
      "dMean.xyio", "dSd.xyio", "dSites.xyio", "dMin.xyio", "dMax.xyio",
      "dTotal.xyio",
      #xy_xi
      "deltaMean.xy_xi", "deltaSum.xy_xi", "deltaRatio.xy_xi",
      "deltaMean_.xy_xi", "deltaSum_.xy_xi", "deltaRatio_.xy_xi",
      "deltaMin.xy_xi", "deltaMinSum.xy_xi", "deltaMinRatio.xy_xi",
      #yi_xy
      "deltaMean.yi_xy", "deltaSum.yi_xy", "deltaRatio.yi_xy",
      "deltaMean_.yi_xy", "deltaSum_.yi_xy", "deltaRatio_.yi_xy",
      "deltaMin.yi_xy", "deltaMinSum.yi_xy", "deltaMinRatio.yi_xy",
      #yi_xi
      "deltaMean.yi_xi", "deltaSum.yi_xi", "deltaRatio.yi_xi",
      "deltaMean_.yi_xi", "deltaSum_.yi_xi", "deltaRatio_.yi_xi",
      "deltaMin.yi_xi", "deltaMinSum.yi_xi", "deltaMinRatio.yi_xi",
      #
      "RND.xyi", "Gmin.xyi", "RNDmin.xyi",
      #
      "RND.xyo", "Gmin.xyo", "RNDmin.xyo")
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
    #xi
    OUT$dMean.xi<-mean(dIUPAC$distIUPAC[x.pos, i.pos], na.rm=TRUE)
    OUT$dMean_.xi<-OUT$dMean.xi - (OUT$dMean.x/2) - (OUT$dMean.i/2)
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
    #yi
    OUT$dMean.yi<-mean(dIUPAC$distIUPAC[y.pos, i.pos], na.rm=TRUE)
    OUT$dMean_.yi<-OUT$dMean.yi - (OUT$dMean.y/2) - (OUT$dMean.i/2)
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
    #io
    OUT$dMean.io<-mean(dIUPAC$distIUPAC[i.pos, o.pos], na.rm=TRUE)
    OUT$dMean_.io<-OUT$dMean.io - (OUT$dMean.i/2) - (OUT$dMean.o/2)
    OUT$dSd.io<-sd(dIUPAC$distIUPAC[i.pos, o.pos], na.rm=TRUE)
    OUT$dSites.io<-mean(dIUPAC$sitesUsed[i.pos, o.pos], na.rm=TRUE)
    OUT$dMin.io<-min(dIUPAC$distIUPAC[i.pos, o.pos], na.rm=TRUE)
    OUT$dMax.io<-max(dIUPAC$distIUPAC[i.pos, o.pos], na.rm=TRUE)
    OUT$dTotal.io<-mean(as.dist(dIUPAC$distIUPAC[c(i.pos, o.pos),
      c(i.pos, o.pos)]), na.rm=TRUE)
    OUT$dSweighted.io<-((length(i.pos)/(length(c(i.pos, o.pos)))) *
      OUT$dMean.i) + ((length(o.pos)/(length(c(i.pos, o.pos)))) * OUT$dMean.o)
    OUT$Fst.io<-(OUT$dTotal.io - OUT$dSweighted.io) / OUT$dTotal.io
    OUT$dRelative.io<-OUT$dMean.io - OUT$dSweighted.io
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
    #xyo
    OUT$dMean.xyo<-mean(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dSd.xyo<-sd(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dSites.xyo<-mean(c(dIUPAC$sitesUsed[x.pos, c(y.pos, o.pos)],
      dIUPAC$sitesUsed[y.pos, o.pos]), na.rm=TRUE)
    OUT$dMin.xyo<-min(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dMax.xyo<-max(c(dIUPAC$distIUPAC[x.pos, c(y.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, o.pos]), na.rm=TRUE)
    OUT$dTotal.xyo<-mean(as.dist(dIUPAC$distIUPAC[c(x.pos, y.pos, o.pos), 
      c(x.pos, y.pos, o.pos)]), na.rm=TRUE)
    #yio
    OUT$dMean.yio<-mean(c(dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dSd.yio<-sd(c(dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dSites.yio<-mean(c(dIUPAC$sitesUsed[y.pos, c(i.pos, o.pos)],
      dIUPAC$sitesUsed[i.pos, o.pos]), na.rm=TRUE)
    OUT$dMin.yio<-min(c(dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dMax.yio<-max(c(dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dTotal.yio<-mean(as.dist(dIUPAC$distIUPAC[c(y.pos, i.pos, o.pos), 
      c(y.pos, i.pos, o.pos)]), na.rm=TRUE)
    #xyio
    OUT$dMean.xyio<-mean(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dSd.xyio<-sd(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dSites.xyio<-mean(c(dIUPAC$sitesUsed[x.pos, c(y.pos, i.pos, o.pos)],
      dIUPAC$sitesUsed[y.pos, c(i.pos, o.pos)],
      dIUPAC$sitesUsed[i.pos, o.pos]), na.rm=TRUE)
    OUT$dMin.xyio<-min(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dMax.xyio<-max(c(dIUPAC$distIUPAC[x.pos, c(y.pos, i.pos, o.pos)],
      dIUPAC$distIUPAC[y.pos, c(i.pos, o.pos)],
      dIUPAC$distIUPAC[i.pos, o.pos]), na.rm=TRUE)
    OUT$dTotal.xyio<-mean(as.dist(dIUPAC$distIUPAC[
      c(x.pos, y.pos, i.pos, o.pos), 
      c(x.pos, y.pos, i.pos, o.pos)]), na.rm=TRUE)
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
    #
    OUT$RND.xyo<-OUT$dMean.xy/((OUT$dMean.xo + OUT$dMean.yo)/2)
    OUT$Gmin.xyo<-OUT$dMin.xy/OUT$dMean.xy
    OUT$RNDmin.xyo<-OUT$dMin.xy/((OUT$dMean.xo + OUT$dMean.yo)/2)
    return(OUT)
}
