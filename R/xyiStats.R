#' @title xyiStats
#' @name xyiStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver and y: donor) with an ingroup
#' population (i: ingroup).
#' In the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined for deltaMean and deltaMin statistics as
#' follows [x:P2 y:P3 i:P1].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined for RND, Gmin and RNDmin
#' statistics as follows [x:P3 y:P2 i:P1].
#' To calcultae RND, Gmin and RNDmin see \code{\link[distIUPAC]{xyoStats}}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [P2 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param y.pos population Y positions [P3 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param i.pos population I positions [P1 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param i.name population I name [default: "i"]
#' @param chr.name chromosome name  [default: "chr"]
#' @param wlen sliding windows length [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites} [default: "bp"]
#' @param dist distance to use [default: "IUPAC"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global or pairwise way (default is to delete in a
#' global way) [default: TRUE]
#' @param threads number of parallel threads [default: 1]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \code{\link[distIUPAC]{rcpp_distIUPAC}}
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @examples
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' AFG.SPRE.CAS.xyiStats<-xyiStats(MySequences, x.pos=AFG.pos,
#' y.pos=SPRE.pos, i.pos=CAS.pos, x.name="AFG", y.name="SPRE", i.name="CAS",
#' threads=2)
#' AFG.SPRE.CAS.xyiStats
#' @export xyiStats
#' @author Kristian K Ullrich
xyiStats<-function(dna, x.pos, y.pos, i.pos, x.name="x", y.name="y",
  i.name="i", chr.name="chr", wlen=25000, wjump=25000, start.by=1, end.by=NULL,
  wtype="bp", dist="IUPAC", global.deletion=TRUE, threads=1, ncores=1,
  pB=FALSE){
    options(scipen=22)
    dna_<-dna[c(x.pos, y.pos, i.pos)]
    x.pos_<-seq(1, length(x.pos))
    y.pos_<-seq(length(x.pos_) + 1, length(x.pos_) + length(y.pos))
    i.pos_<-seq(length(x.pos_) + length(y.pos_) + 1,
      length(x.pos_) + length(y.pos_) + length(i.pos))
    #OUT<-distIUPACsw(dna_, FUN=function(x) {
    #      dist.xyiStats(x, x.pos=x.pos_, y.pos=y.pos_, i.pos=i.pos_,
    #      x.name=x.name, y.name=y.name, i.name=i.name)
    #  }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
    #  end.by=end.by, wtype=wtype, dist=dist, global.deletion=global.deletion,
    #  threads=threads, ncores=ncores, pB=pB)
    OUT<-tmpSEQsw(dna_, FUN=function(x) {
          dist.xyiStats(x, x.pos=x.pos_, y.pos=y.pos_, i.pos=i.pos_,
          x.name=x.name, y.name=y.name, i.name=i.name, dist=dist, ncores=ncores)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
      end.by=end.by, wtype=wtype, global.deletion=global.deletion,
      threads=threads, pB=pB)
    return(OUT)
}
