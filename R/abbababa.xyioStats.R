#' @title abbababa.xyioStats
#' @name abbababa.xyioStats
#' @description This function calculates ABBA-BABA statistics
#' comparing two populations (x: receiver; y: donor) with an ingroup
#' population (i: ingroup) and an outgroup population (o: outgroup).
#' In a four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined as follows [x:P2 y:P3 i:P1 o:P4].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined as follows [x:P3 y:P2 i:P1 o:P4].
#' @import foreach
#' @importFrom stats as.dist sd setNames
#' @param tmpSEQ \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [mandatory]
#' @param y.pos population Y positions [mandatory]
#' @param i.pos population I positions [mandatory]
#' @param o.pos population O positions [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param i.name population I name [default: "i"]
#' @param o.name population O name [default: "o"]
#' @param x.freq minimal frequency for population X to keep site [default: 1.0]
#' @param y.freq minimal frequency for population Y to keep site [default: 1.0]
#' @param i.freq minimal frequency for population I to keep site [default: 1.0]
#' @param o.freq minimal frequency for population O to keep site [default: 1.0]
#' @seealso \code{\link[distIUPAC]{xyioStats}},
#' \code{\link[distIUPAC]{diploid2haploid}}
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
  x.name="x", y.name="y", i.name="i", o.name="o"){
    #forked from https://github.com/simonhmartin/genomics_general
    calc.f4<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        ((1 - p1dFreq) * p2dFreq * p3dFreq * (1 - p4dFreq)) - 
        (p1dFreq * (1 - p2dFreq) * p3dFreq * (1 - p4dFreq))
    }
    calc.f4_c<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        calc.f4(p1dFreq, p2dFreq, p3dFreq, p4dFreq) +
        calc.f4((1 - p1dFreq), (1 - p2dFreq), (1 - p3dFreq), (1 - p4dFreq))
    }
    calc.ABAAsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (1 - p1dFreq) * p2dFreq * (1 - p3dFreq) * (1 - p4dFreq)
    }
    calc.BAAAsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        p1dFreq * (1 - p2dFreq) * (1 - p3dFreq) * (1 - p4dFreq)
    }
    calc.ABBAsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (1 - p1dFreq) * p2dFreq * p3dFreq * (1 - p4dFreq)
    }
    calc.BABAsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        p1dFreq * (1 - p2dFreq) * p3dFreq * (1 - p4dFreq)
    }
    calc.ABAA_BABBsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        ((1 - p1dFreq) * p2dFreq * (1 - p3dFreq) * (1 - p4dFreq)) +
        (p1dFreq * (1 - p2dFreq) * p3dFreq * p4dFreq)
    }
    calc.BAAA_ABBBsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (p1dFreq * (1 - p2dFreq) * (1 - p3dFreq) * (1 - p4dFreq)) +
        ((1 - p1dFreq) * p2dFreq * p3dFreq * p4dFreq)
    }
    calc.ABBA_BAABsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        ((1 - p1dFreq) * p2dFreq * p3dFreq * (1 - p4dFreq)) +
        (p1dFreq * (1 - p2dFreq) * (1 - p3dFreq) * p4dFreq)
    }
    calc.BABA_ABABsum<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (p1dFreq * (1 - p2dFreq) * p3dFreq * (1 - p4dFreq)) +
        ((1 - p1dFreq) * p2dFreq * (1 - p3dFreq) * p4dFreq)
    }
    calc.maxABBAsumHom<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (1 - p1dFreq) * p3dFreq * p3dFreq * (1 - p4dFreq)
    }
    calc.maxBABAsumHom<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        p1dFreq * (1 - p3dFreq) * p3dFreq * (1 - p4dFreq)
    }
    calc.maxABBAsumD<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        OUT<-NULL
        foreach(i=seq_along(p1dFreq)) %do% {
            if(p3dFreq[i] >= p2dFreq[i]){
                OUT<-c(OUT,
                  (1 - p1dFreq[i]) * p3dFreq[i] * p3dFreq[i] * (1 - p4dFreq[i]))
            }
            else{
                OUT<-c(OUT,
                  (1 - p1dFreq[i]) * p2dFreq[i] * p2dFreq[i] * (1 - p4dFreq[i]))
            }
        }
        return(OUT)
    }
    calc.maxBABAsumD<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        OUT<-NULL
        foreach(i=seq_along(p1dFreq)) %do% {
            if(p3dFreq[i] >= p2dFreq[i]){
                OUT<-c(OUT,
                  p1dFreq[i] * (1 - p3dFreq[i]) * p3dFreq[i] * (1 - p4dFreq[i]))
            }
            else{
                OUT<-c(OUT,
                  p1dFreq[i] * (1 - p2dFreq[i]) * p2dFreq[i] * (1 - p4dFreq[i]))
            }
        }
        return(OUT)
    }
    calc.maxABBAsumG<-function(p1dFreq, p3adFreq, p3bdFreq, p4dFreq){
        (1 - p1dFreq) * p3adFreq * p3bdFreq * (1 - p4dFreq)
    }
    calc.maxBABAsumG<-function(p1dFreq, p3adFreq, p3bdFreq, p4dFreq){
        p1dFreq * (1 - p3adFreq) * p3bdFreq * (1 - p4dFreq)
    }
    calc.ABBAsumG<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        (1 - p1dFreq) * p2dFreq * p3dFreq * (1 - p4dFreq)
    }
    calc.BABAsumG<-function(p1dFreq, p2dFreq, p3dFreq, p4dFreq){
        p1dFreq * (1 - p2dFreq) * p3dFreq * (1 - p4dFreq)
    }
    options(scipen=22)
    calc.real<-FALSE
    if(length(y.pos)>1){
        calc.real<-TRUE
        ya.pos<-y.pos[1:floor(length(y.pos)/2)]
        yb.pos<-y.pos[(floor(length(y.pos)/2)+1):length(y.pos)]
    }
    XNAME<-x.name
    YNAME<-y.name
    INAME<-i.name
    ONAME<-o.name
    #freq
    freq.x<-NA
    freq.y<-NA
    freq.i<-NA
    freq.o<-NA
    freq.x.rem<-NA
    freq.y.rem<-NA
    freq.i.rem<-NA
    freq.o.rem<-NA
    freq.all.rem<-NA
    #
    ABBAsum<-0
    BABAsum<-0
    ABBAsumG<-0
    BABAsumG<-0
    maxABBAsumG<-0
    maxBABAsumG<-0
    maxABBAsumHom<-0
    maxBABAsumHom<-0
    maxABBAsumD<-0
    maxBABAsumD<-0
    #
    remcount<-0
    remcount.x<-0
    remcount.y<-0
    remcount.i<-0
    remcount.o<-0
    monocount<-0
    bicount<-0
    tricount<-0
    tetracount<-0
    #
    D<-NA
    fG<-NA
    fhom<-NA
    fd<-NA
    OUT<-list(XNAME, YNAME, INAME, ONAME,
      ABBAsum, BABAsum, ABBAsumG, BABAsumG, maxABBAsumG, maxBABAsumG,
      maxABBAsumHom, maxBABAsumHom, maxABBAsumD, maxBABAsumD,
      monocount, bicount, tricount, tetracount,
      remcount, remcount.x, remcount.y, remcount.i, remcount.o,
      D, fG, fhom, fd)
    names(OUT)<-c("XNAME", "YNAME", "INAME", "ONAME",
      "ABBAsum", "BABAsum", "ABBAsumG", "BABAsumG",
      "maxABBAsumG", "maxBABAsumG", "maxABBAsumHom", "maxBABAsumHom",
      "maxABBAsumD", "maxBABAsumD",
      "monocount", "bicount", "tricount", "tetracount",
      "remcount", "remcount.x", "remcount.y", "remcount.i", "remcount.o",
      "D", "fG","fhom", "fd")
    #
    cM<-consensusMatrix(iupac2diploid(tmpSEQ))
    x.cM<-consensusMatrix(iupac2diploid(tmpSEQ[x.pos]))
    y.cM<-consensusMatrix(iupac2diploid(tmpSEQ[y.pos]))
    if(calc.real){
        ya.cM<-consensusMatrix(iupac2diploid(tmpSEQ[ya.pos]))
        yb.cM<-consensusMatrix(iupac2diploid(tmpSEQ[yb.pos]))
    }
    i.cM<-consensusMatrix(iupac2diploid(tmpSEQ[i.pos]))
    o.cM<-consensusMatrix(iupac2diploid(tmpSEQ[o.pos]))
    #get frequencies
    freq.x<-colSums(x.cM[1:4, ])/(2*length(x.pos))
    freq.y<-colSums(y.cM[1:4, ])/(2*length(y.pos))
    freq.i<-colSums(i.cM[1:4, ])/(2*length(i.pos))
    freq.o<-colSums(o.cM[1:4, ])/(2*length(o.pos))
    #get idx of sites to be reomved based on minimal frequency per population
    freq.x.rem<-which(freq.x<x.freq)
    freq.y.rem<-which(freq.y<y.freq)
    freq.i.rem<-which(freq.i<i.freq)
    freq.o.rem<-which(freq.o<o.freq)
    freq.all.rem<-unique(c(freq.x.rem, freq.y.rem, freq.i.rem, freq.o.rem))
    OUT$remcount<-length(freq.all.rem)
    OUT$remcount.x<-length(freq.x.rem)
    OUT$remcount.y<-length(freq.y.rem)
    OUT$remcount.i<-length(freq.i.rem)
    OUT$remcount.o<-length(freq.o.rem)
    #reduce to sites that are kept
    if(length(freq.all.rem)!=0){
        cM.red<-cM[, -freq.all.rem, drop=FALSE]
        x.cM.red<-x.cM[, -freq.all.rem, drop=FALSE]
        y.cM.red<-y.cM[, -freq.all.rem, drop=FALSE]
        i.cM.red<-i.cM[, -freq.all.rem, drop=FALSE]
        o.cM.red<-o.cM[, -freq.all.rem, drop=FALSE]
        if(calc.real){
            ya.cM.red<-ya.cM[, -freq.all.rem, drop=FALSE]
            yb.cM.red<-yb.cM[, -freq.all.rem, drop=FALSE]
        }
    }
    if(length(freq.all.rem)==0){
        cM.red<-cM
        x.cM.red<-x.cM
        y.cM.red<-y.cM
        i.cM.red<-i.cM
        o.cM.red<-o.cM
        if(calc.real){
            ya.cM.red<-ya.cM
            yb.cM.red<-yb.cM
        }
    }
    #get allele counts
    alleles<-apply(cM.red[1:4,], 2, function(x) length(which(x>0)))
    OUT$monocount<-length(which(alleles==1))
    OUT$bicount<-length(which(alleles==2))
    OUT$tricount<-length(which(alleles==3))
    OUT$tetracount<-length(which(alleles==4))
    #consider only bicount sites
    cM.red.bi<-cM.red[, alleles==2, drop=FALSE]
    x.cM.red.bi<-x.cM.red[, alleles==2, drop=FALSE]
    y.cM.red.bi<-y.cM.red[, alleles==2, drop=FALSE]
    i.cM.red.bi<-i.cM.red[, alleles==2, drop=FALSE]
    o.cM.red.bi<-o.cM.red[, alleles==2, drop=FALSE]
    if(calc.real){
        ya.cM.red.bi<-ya.cM.red[, alleles==2, drop=FALSE]
        yb.cM.red.bi<-yb.cM.red[, alleles==2, drop=FALSE]
    }
    #get derived and ancestral state
    #if the outgroup is fixed, then that is the ancestral state
    #- otherwise the anc state is the most common allele overall
    #get fixed outgroup positions
    o.alleles<-apply(o.cM.red.bi, 2, function(x) length(which(x>0)))
    o.fixed<-which(o.alleles==1)
    if(length(o.fixed)>0){
        ancestral.fixed<-apply(o.cM.red.bi[1:4, o.fixed, drop=FALSE], 2,
          function(x) which(x==max(x)))
        derived.fixed<-apply(rbind(cM.red.bi[1:4, o.fixed, drop=FALSE],
          ancestral.fixed), 2, function(x) which(x[1:4]>0)
          [which(x[1:4]>0)!=x[5]])
        x.derived.fixed.freq<-(x.cM.red.bi[
          cbind(derived.fixed, o.fixed)])/
          (colSums(x.cM.red.bi[1:4, o.fixed]))
        y.derived.fixed.freq<-(y.cM.red.bi[
          cbind(derived.fixed, o.fixed)])/
          (colSums(y.cM.red.bi[1:4, o.fixed]))
        i.derived.fixed.freq<-(i.cM.red.bi[
          cbind(derived.fixed, o.fixed)])/
          (colSums(i.cM.red.bi[1:4, o.fixed]))
        o.derived.fixed.freq<-(o.cM.red.bi[
          cbind(derived.fixed, o.fixed)])/
          (colSums(o.cM.red.bi[1:4, o.fixed]))
        if(calc.real){
            ya.derived.fixed.freq<-(ya.cM.red.bi[
              cbind(derived.fixed, o.fixed)])/
              (colSums(ya.cM.red.bi[1:4, o.fixed]))
            yb.derived.fixed.freq<-(yb.cM.red.bi[
              cbind(derived.fixed, o.fixed)])/
              (colSums(yb.cM.red.bi[1:4, o.fixed]))
        }
        ABBAsum<-sum(ABBAsum, calc.ABBAsum(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        BABAsum<-sum(BABAsum, calc.BABAsum(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        maxABBAsumHom<-sum(maxABBAsumHom, calc.maxABBAsumHom(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        maxBABAsumHom<-sum(maxBABAsumHom, calc.maxBABAsumHom(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        maxABBAsumD<-sum(maxABBAsumD, calc.maxABBAsumD(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        maxBABAsumD<-sum(maxBABAsumD, calc.maxBABAsumD(
          i.derived.fixed.freq,
          x.derived.fixed.freq,
          y.derived.fixed.freq,
          o.derived.fixed.freq))
        if(calc.real){
            maxABBAsumG<-sum(maxABBAsumG, calc.maxABBAsumG(
              i.derived.fixed.freq,
              ya.derived.fixed.freq,
              yb.derived.fixed.freq,
              o.derived.fixed.freq))
            maxBABAsumG<-sum(maxBABAsumG, calc.maxBABAsumG(
              i.derived.fixed.freq,
              ya.derived.fixed.freq,
              yb.derived.fixed.freq,
              o.derived.fixed.freq))
            ABBAsumG<-sum(ABBAsumG, calc.ABBAsumG(
              i.derived.fixed.freq,
              x.derived.fixed.freq,
              y.derived.fixed.freq,
              o.derived.fixed.freq))
            BABAsumG<-sum(BABAsumG, calc.BABAsumG(
              i.derived.fixed.freq,
              x.derived.fixed.freq,
              y.derived.fixed.freq,
              o.derived.fixed.freq))
        }
    }
    #get most common state for unfixed outgroup positions
    o.unfixed<-which(o.alleles!=1)
    if(length(o.unfixed)>0){
        ancestral.unfixed<-apply(cM.red.bi[1:4, o.unfixed, drop=FALSE], 2,
          function(x) which(x==max(x)))
        derived.unfixed<-apply(cM.red.bi[1:4, o.unfixed, drop=FALSE], 2,
          function(x) which(x==min(x[x>0])))
        #get derived unfixed frequencies per population
        x.derived.unfixed.freq<-(x.cM.red.bi[
          cbind(derived.unfixed, o.unfixed)])/
          (colSums(x.cM.red.bi[1:4, o.unfixed]))
        y.derived.unfixed.freq<-(y.cM.red.bi[
          cbind(derived.unfixed, o.unfixed)])/
          (colSums(y.cM.red.bi[1:4, o.unfixed]))
        i.derived.unfixed.freq<-(i.cM.red.bi[
          cbind(derived.unfixed, o.unfixed)])/
          (colSums(i.cM.red.bi[1:4, o.unfixed]))
        o.derived.unfixed.freq<-(o.cM.red.bi[
          cbind(derived.unfixed, o.unfixed)])/
          (colSums(o.cM.red.bi[1:4, o.unfixed]))
        if(calc.real){
            ya.derived.unfixed.freq<-(ya.cM.red.bi[
              cbind(derived.unfixed, o.unfixed)])/
              (colSums(ya.cM.red.bi[1:4, o.unfixed]))
            yb.derived.unfixed.freq<-(yb.cM.red.bi[
              cbind(derived.unfixed, o.unfixed)])/
              (colSums(yb.cM.red.bi[1:4, o.unfixed]))
        }
        ABBAsum<-sum(ABBAsum, calc.ABBAsum(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        BABAsum<-sum(BABAsum, calc.BABAsum(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        maxABBAsumHom<-sum(maxABBAsumHom, calc.maxABBAsumHom(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        maxBABAsumHom<-sum(maxBABAsumHom, calc.maxBABAsumHom(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        maxABBAsumD<-sum(maxABBAsumD, calc.maxABBAsumD(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        maxBABAsumD<-sum(maxBABAsumD, calc.maxBABAsumD(
          i.derived.unfixed.freq,
          x.derived.unfixed.freq,
          y.derived.unfixed.freq,
          o.derived.unfixed.freq))
        if(calc.real){
            maxABBAsumG<-sum(maxABBAsumG, calc.maxABBAsumG(
              i.derived.unfixed.freq,
              ya.derived.unfixed.freq,
              yb.derived.unfixed.freq,
              o.derived.unfixed.freq))
            maxBABAsumG<-sum(maxBABAsumG, calc.maxBABAsumG(
              i.derived.unfixed.freq,
              ya.derived.unfixed.freq,
              yb.derived.unfixed.freq,
              o.derived.unfixed.freq))
            ABBAsumG<-sum(ABBAsumG, calc.ABBAsumG(
              i.derived.unfixed.freq,
              x.derived.unfixed.freq,
              y.derived.unfixed.freq,
              o.derived.unfixed.freq))
            BABAsumG<-sum(BABAsumG, calc.BABAsumG(
              i.derived.unfixed.freq,
              x.derived.unfixed.freq,
              y.derived.unfixed.freq,
              o.derived.unfixed.freq))
        }
    }
    OUT$ABBAsum<-ABBAsum
    OUT$BABAsum<-BABAsum
    OUT$ABBAsumG<-ABBAsumG
    OUT$BABAsumG<-BABAsumG
    OUT$maxABBAsumG<-maxABBAsumG
    OUT$maxBABAsumG<-maxBABAsumG
    OUT$maxABBAsumHom<-maxABBAsumHom
    OUT$maxBABAsumHom<-maxBABAsumHom
    OUT$maxABBAsumD<-maxABBAsumD
    OUT$maxBABAsumD<-maxBABAsumD
    OUT$D<-(ABBAsum - BABAsum) / (ABBAsum + BABAsum)
    OUT$fG<-(ABBAsumG - BABAsumG) / (maxABBAsumG - maxBABAsumG)
    OUT$fhom<-(ABBAsum - BABAsum) / (maxABBAsumHom - maxBABAsumHom)
    OUT$fd<-(ABBAsum - BABAsum) / (maxABBAsumD - maxBABAsumD)
    return(OUT)
}
