#' @title bam2consensus
#' @name bam2consensus
#' @description This function calculates \code{distIUPAC} based distances
#' within one population (x: receiver; x: donor) along sliding window.
#' @param pop.list population list
#' @param pop.names population names [default: NULL]
#' @param chr.name chromosome name [default: NULL]
#' @param wlen sliding windows length [default: 100000]
#' @param wjump sliding windows jump [default: 100000]
#' @param start.by optional start position [default: NULL]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites} [default: "bp"] see \link[distIUPAC]{biSites} and
#' \link[distIUPAC]{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global [TRUE] or pairwise way [FALSE]
#' @param threads number of parallel threads to process sliding windows
#' @param ncores number of parallel cores to process pairwise
#' distance calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC}
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' ##multiple threads to process 
#' CAS.xStats<-xStats(MySequences, x.pos=CAS.pos, x.name="CAS", threads=2)
#' CAS.xStats
#' @export bam2consensus
#' @author Kristian K Ullrich
bam2consensus<-function(dna, x.pos=NULL, x.name="x", chr.name="chr",
  wlen=25000, wjump=25000, start.by=1, end.by=NULL, wtype="bp",
  dist="IUPAC", global.deletion=TRUE, threads=1, ncores=1, pB=TRUE){
    options(scipen=22)
    if(is.null(x.pos)){
        x.pos<-seq(1,length(dna))
    }
    dna_<-dna[x.pos]
    x.pos_<-seq(1,length(x.pos))
    OUT<-distIUPACsw(dna_,
      FUN=function(x) {
          dist.xStats(x, x.pos=x.pos_, x.name=x.name)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump,
      start.by=start.by, end.by=end.by, wtype=wtype, dist=dist,
      global.deletion=global.deletion, threads=threads, ncores=ncores, pB=pB)
    return(OUT)
}
