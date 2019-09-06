#' @title dist.xStats
#' @name dist.xStats
#' @description This function calculates \code{distIUPAC} based distances
#' within one population (x: receiver; x: donor).
#' @importFrom stats as.dist sd setNames
#' @param tmpSEQ \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [default: NULL]
#' @param x.name population X name [default: "x"]
#' @param dist distance to use [default: IUPAC] or choose one model as in
#' \link[ape]{dist.dna} [default: "IUPAC"]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC} [default: 1]
#' @seealso \code{\link[distIUPAC]{xStats}},
#' \code{\link[distIUPAC]{distIUPAC}},
#' \code{\link[distIUPAC]{rcpp_distIUPAC}}
#' @examples
#' ##Use the 'xStats' function to handle population assignment automatically
#' 
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' 
#' ##pairwise deletion
#' CAS.pairwiseDeletion.dist<-distIUPAC(as.character(MySequences[CAS.pos]))
#' CAS.pairwiseDeletion.xStats<-dist.xStats(
#'   CAS.pairwiseDeletion.dist, x.name="CAS")
#' 
#' ##global deletion
#' CAS.globalDeletion.dist<-distIUPAC(
#'   as.character(globalDeletion(MySequences[CAS.pos])))
#' CAS.globalDeletion.xStats<-dist.xStats(
#'   CAS.globalDeletion.dist, x.name="CAS")
#' 
#' ##compare results
#' rbind(CAS.pairwiseDeletion.xStats, CAS.globalDeletion.xStats)
#' @export dist.xStats
#' @author Kristian K Ullrich
#dist.xStats<-function(dIUPAC, x.pos=NULL, x.name="x"){
dist.xStats<-function(tmpSEQ, x.pos=NULL,
  x.name="x", dist="IUPAC", ncores=1){
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
    if(is.null(x.pos)){
        x.pos<-seq(from=1, to=ncol(dIUPAC$distIUPAC))
    }
    XNAME<-x.name
    dMean.x<-NA
    dSd.x<-NA
    dSites.x<-NA
    dMin.x<-NA
    dMax.x<-NA
    OUT<-list(XNAME, dMean.x, dSd.x, dSites.x, dMin.x, dMax.x)
    names(OUT)<-c("XNAME", "dMean.x", "dSd.x", "dSites.x", "dMin.x", "dMax.x")
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
    return(OUT)
}
