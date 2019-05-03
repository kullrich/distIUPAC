#' @title xStats
#' @name xStats
#' @description This function calculates \code{distIUPAC} based distances within one population
#' (x: receiver; x: donor).
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param start.by optional start position
#' @param end.by optional end position
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites with missing data in a global or pairwise way (default is to delete in a global way)
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param chr.name chromosome name
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' CAS.xStats<-xStats(MySequences, x.pos = CAS.pos, x.name = "CAS", threads = 2)
#' CAS.xStats
#' @export xStats
#' @author Kristian K Ullrich
xStats<-function(dna, x.pos=NULL, x.name="x", chr.name="chr",
 wlen=25000, wjump=25000, start.by=NULL, end.by=NULL, wtype="bp",
 dist="IUPAC", global.deletion=TRUE, threads=1, pB=TRUE){
  options(scipen=22)
  if(is.null(x.pos)){
    x.pos<-seq(1,length(dna))
  }
  dna_<-dna[x.pos]
  x.pos_<-seq(1,length(x.pos))
  OUT<-distIUPACsw(dna_, FUN=function(x) {dist.xStats(x, x.pos=x.pos_, x.name=x.name)},
   chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by, end.by=end.by, wtype=wtype,
   dist=dist, global.deletion=global.deletion, threads=threads, pB=pB)
  return(OUT)
}
