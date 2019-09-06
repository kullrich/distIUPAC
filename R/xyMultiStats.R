#' @title xyMultiStats
#' @name xyMultiStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing all possible pairwise population combinations (x: receiver and
#' y: donor).
#' @importFrom utils combn
#' @importFrom stats as.dist sd setNames
#' @param dna \code{DNAStringSet} [mandatory]
#' @param pop.list population list with individual positions [mandatory]
#' @param chr.name chromosome name [default: "chr"]
#' @param wlen sliding windows length [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [deafult: 1]
#' @param end.by optional end position [deafult: NULL]
#' @param wtype sliding windows type to use \code{bp},
#' \code{biSites} or \code{triSites} [default: "bp"]
#' @param dist distance to use [deafult: "IUPAC"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global or pairwise way (default is to delete in a
#' global way) [default: TRUE]
#' @param threads number of parallel threads [default: 1]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC}
#' @param pB specifies if progress should be shown as a progress bar
#' [default: TRUE]
#' @examples
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' pop.list<-setNames(list(CAS.pos, AFG.pos, SPRE.pos), c("CAS", "AFG",
#' "SPRE"))
#' ##sliding windows based on base-pair length
#' CAS.AFG.SPRE.xyMultiStats<-xyMultiStats(MySequences, pop.list=pop.list,
#' threads = 2)
#' CAS.AFG.SPRE.xyMultiStats
#' #sliding windows based on biSites
#' CAS.AFG.SPRE.xyMultiStats<-xyMultiStats(MySequences, list.pos = pop.list,
#' wtype = "biSites", wlen = 50, threads = 2)
#' CAS.AFG.SPRE.xyMultiStats
#' @export xyMultiStats
#' @author Kristian K Ullrich
xyMultiStats<-function(dna, pop.list, chr.name="chr", wlen=25000, wjump=25000,
  start.by=1, end.by=NULL, wtype="bp", dist="IUPAC", global.deletion=TRUE,
  threads=1, ncores=1, pB=TRUE){
    options(scipen=22)
    if(!is.list(pop.list)){
        stop("pop.list needs to be a list of positions")
    }
    if(is.null(names(pop.list))){
        names(pop.list)<-seq(from=1,
        to=length(pop.list))
    }
    pop.comb<-combn(length(pop.list), 2)
    COMBOUT<-vector("list", length=dim(pop.comb)[2])
    names(COMBOUT)<-apply(pop.comb, 2, function(x) {
        paste0(x, collapse="_")
    })
    for(c.idx in seq(from=1, to=dim(pop.comb)[2])){
        x.pos<-unlist(pop.list[pop.comb[1, c.idx]])
        y.pos<-unlist(pop.list[pop.comb[2, c.idx]])
        x.name<-names(pop.list[pop.comb[1, c.idx]])
        y.name<-names(pop.list[pop.comb[2, c.idx]])
        dna_<-dna[c(x.pos,y.pos)]
        x.pos_<-seq(from=1, to=length(x.pos))
        y.pos_<-seq(from=length(x.pos_) + 1, to=length(x.pos_) + length(y.pos))
        COMBOUT[c.idx]<-list(xyStats(dna=dna_, x.pos=x.pos_, y.pos=y.pos_,
          x.name=x.name, y.name=y.name, chr.name=chr.name, wlen=wlen,
          wjump=wjump, start.by=start.by, end.by=end.by, wtype=wtype,
          dist=dist, global.deletion=global.deletion, threads=threads,
          ncores=ncores, pB=pB))
    }
    return(COMBOUT)
}
