#' @title xyioStats
#' @name xyioStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver and y: donor) with an ingroup
#' population (i: ingroup) and an outgroup population (o: outgroup).
#' In the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined for deltaMean and deltaMin statistics as
#' follows [x:P2 y:P3 i:P1 o:P4].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined for RND, Gmin and RNDmin
#' statistics as follows [x:P3 y:P2 i:P1 o:P4]. Optional, ABBA-BABA statisctis
#' can be calculated concordantly on the four-taxon scenario
#' (((P1:i,P2:x),P3:y),O:o).
#' @param dna \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [P2 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param y.pos population Y positions [P3 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param i.pos population I positions [P1 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param o.pos population I positions [P4 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param i.name population I name [default: "i"]
#' @param o.name population O name [default: "o"]
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
#' @param do.ABBA specifies if ABBA-BABA statistics should be calculated
#' [default: FALSE]
#' @param x.freq [default: 1.0]
#' @param y.freq [default: 1.0]
#' @param i.freq [default: 1.0]
#' @param o.freq [default: 1.0]
#' @examples
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' APO.pos<-1
#' AFG.SPRE.CAS.APO.xyioStats<-xyioStats(MySequences, x.pos=AFG.pos,
#' y.pos=SPRE.pos, i.pos=CAS.pos, o.pos=APO.pos,
#' x.name="AFG", y.name="SPRE", i.name="CAS", o.name="APO",
#' threads=2)
#' AFG.SPRE.CAS.APO.xyioStats
#' @export xyioStats
#' @author Kristian K Ullrich
xyioStats<-function(dna, x.pos, y.pos, i.pos, o.pos, x.name="x", y.name="y",
  i.name="i", o.name="o", chr.name="chr", wlen=25000, wjump=25000, start.by=1,
  end.by=NULL, wtype="bp", dist="IUPAC", global.deletion=TRUE, threads=1,
  ncores=1, pB=FALSE, do.ABBA=FALSE, x.frq=1.0, y.freq=1.0, i.freq=1.0,
  o.freq=1.0){
    options(scipen=22)
    dna_<-dna[c(x.pos, y.pos, i.pos, o.pos)]
    x.pos_<-seq(1, length(x.pos))
    y.pos_<-seq(length(x.pos_) + 1, length(x.pos_) + length(y.pos))
    i.pos_<-seq(length(x.pos_) + length(y.pos_) + 1,
      length(x.pos_) + length(y.pos_) + length(i.pos))
    o.pos_<-seq(length(x.pos_) + length(y.pos_) + length(i.pos_) + 1,
      length(x.pos_) + length(y.pos_) + length(i.pos_) + length(o.pos))
    #OUT<-distIUPACsw(dna_, FUN=function(x) {
    #      dist.xyioStats(x, x.pos=x.pos_, y.pos=y.pos_, i.pos=i.pos_,
    #      o.pos=o.pos_, x.name=x.name, y.name=y.name, i.name=i.name,
    #      o.name=o.name)
    #  }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
    #  end.by=end.by, wtype=wtype, dist=dist, global.deletion=global.deletion,
    #  threads=threads, ncores=ncores, pB=pB)
    OUT<-tmpSEQsw(dna_, FUN=function(x) {
          dist.xyioStats(x, x.pos=x.pos_, y.pos=y.pos_, i.pos=i.pos_,
          o.pos=o.pos_, x.name=x.name, y.name=y.name, i.name=i.name,
          o.name=o.name, dist=dist, ncores=ncores)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
      end.by=end.by, wtype=wtype, global.deletion=global.deletion,
      threads=threads, pB=pB)
    if(do.ABBA){
        ABBA<-tmpSEQsw(dna_, FUN=function(x) {
              abbababa.xyioStats(x, x.pos=x.pos_, y.pos=y.pos_, i.pos=i.pos_,
              o.pos=o.pos_, x.name=x.name, y.name=y.name, i.name=i.name,
              o.name=o.name, x.freq=x.freq, y.freq=y.freq, i.freq=i.freq,
              o.freq=o.freq, ncores=ncores)
          }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
          end.by=end.by, wtype=wtype, global.deletion=global.deletion,
          threads=threads, pB=pB)
        return(setNames(list(OUT,ABBA),c("xyioStats","xyioABBA")))
    }
    return(OUT)
}
