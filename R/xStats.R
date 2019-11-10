#' @title xStats
#' @name xStats
#' @description This function calculates \code{distIUPAC} based distances
#' within one population (x: receiver; x: donor) along sliding window.
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions [default: NULL]
#' take all sequences from \code{DNAStringSet}
#' @param x.name population X name [default: "x"]
#' @param chr.name chromosome name [default: "chr"]
#' @param wlen sliding windows length  [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites}  [default: "bp"]
#' @param dist distance to use [default: IUPAC] or choose one model as in
#' \link[ape]{dist.dna} [default: "IUPAC"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global [default: TRUE] or pairwise way [FALSE]
#' @param threads number of parallel threads [default: 1]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @references Paradis, E., & Schliep, K. (2018). ape 5.0: an environment for
#' modern phylogenetics and evolutionary analyses in R. \emph{Bioinformatics},
#' \bold{35(3)}, 526-528.
#' @seealso \code{\link[distIUPAC]{dist.xStats}},
#' \code{\link[ape]{dist.dna}}
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' ##multiple threads to process 
#' CAS.xStats<-xStats(MySequences, x.pos=CAS.pos, x.name="CAS", threads=2)
#' CAS.xStats
#' @export xStats
#' @author Kristian K Ullrich
xStats<-function(dna, x.pos=NULL, x.name="x", chr.name="chr",
  wlen=25000, wjump=25000, start.by=1, end.by=NULL, wtype="bp",
  dist="IUPAC", global.deletion=TRUE, threads=1, ncores=1, pB=FALSE){
    options(scipen=22)
    if(is.null(x.pos)){
        x.pos<-seq(1,length(dna))
    }
    dna_<-dna[x.pos]
    x.pos_<-seq(1,length(x.pos))
    #OUT<-distIUPACsw(dna_,
    #  FUN=function(x) {
    #      dist.xStats(x, x.pos=x.pos_, x.name=x.name)
    #  }, chr.name=chr.name, wlen=wlen, wjump=wjump,
    #  start.by=start.by, end.by=end.by, wtype=wtype, dist=dist,
    #  global.deletion=global.deletion, threads=threads, ncores=ncores, pB=pB)
    OUT<-tmpSEQsw(dna_,
      FUN=function(x) {
          dist.xStats(x, x.pos=x.pos_, x.name=x.name, dist=dist, ncores=ncores)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump,
      start.by=start.by, end.by=end.by, wtype=wtype,
      global.deletion=global.deletion, threads=threads, pB=pB)
    return(OUT)
}
