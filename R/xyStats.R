#' @title xyStats
#' @name xyStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver and y: donor).
#' @param dna \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions [mandatory]
#' @param y.pos population Y positions [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param chr.name chromosome name [default: "chr"]
#' @param wlen sliding windows length [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites} [default: "bp"]
#' @param dist distance to use, choose one model as in
#' \link[ape]{dist.dna} or [default: "IUPAC"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global or pairwise way (default is to delete in a
#' global way) [default: TRUE]
#' @param threads number of parallel threads [default: 1]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @references Slatkin, M. (1991). Inbreeding coefficients and coalescence 
#' times. \emph{Genetics Research}, \bold{58(2)}, 167-175.
#' 
#' Beerli, P. (1998). Structured Populations. \emph{Advances in molecular
#' ecology}, \bold{306}, 39.
#' 
#' Reich, D., Thangaraj, K., Patterson, N., Price, A. L., & Singh,
#' L. (2009). Reconstructing Indian population history. \emph{Nature},
#' \bold{461(7263)}, 489.
#' 
#' Patterson, N., Moorjani, P., Luo, Y., Mallick, S., Rohland, N., Zhan, Y.,
#' ... & Reich, D. (2012). Ancient admixture in human history. \emph{Genetics},
#' \bold{192(3)}, 1065-1093.
#' 
#' Peter, B. M. (2016). Admixture, population structure, and F-statistics.
#' \emph{Genetics}, \bold{202(4)}, 1485-1501.
#' 
#' Paradis, E., & Schliep, K. (2018). ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R. \emph{Bioinformatics},
#' \bold{35(3)}, 526-528.
#' @seealso \code{\link[distIUPAC]{dist.xyStats}},
#' \code{\link[ape]{dist.dna}}
#' @examples
#' ##load sequence data
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' ##sliding windows based on base-pair length
#' CAS.AFG.xyStats<-xyStats(MySequences, x.pos=CAS.pos, y.pos=AFG.pos,
#' x.name="CAS", y.name="AFG", threads=2)
#' CAS.AFG.xyStats
#' ##sliding windows based on biSites
#' CAS.AFG.xyStats<-xyStats(MySequences, x.pos=CAS.pos, y.pos=AFG.pos,
#' wlen=50, wtype="biSites", x.name="CAS", y.name="AFG", threads=2)
#' CAS.AFG.xyStats
#' @export xyStats
#' @author Kristian K Ullrich
xyStats<-function(dna, x.pos, y.pos, x.name="x", y.name="y", chr.name="chr",
  wlen=25000, wjump=25000, start.by=1, end.by=NULL, wtype="bp", dist="IUPAC",
  global.deletion=TRUE, threads=1, ncores=1, pB=FALSE){
    options(scipen=22)
    dna_<-dna[c(x.pos,y.pos)]
    x.pos_<-seq(1,length(x.pos))
    y.pos_<-seq(length(x.pos_)+1,length(x.pos_)+length(y.pos))
    #OUT<-distIUPACsw(dna_, FUN=function(x) {
    #      dist.xyStats(x, x.pos=x.pos_, y.pos=y.pos_, x.name=x.name,
    #      y.name=y.name)
    #  }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
    #  end.by=end.by, wtype=wtype, dist=dist, global.deletion=global.deletion,
    #  threads=threads, ncores=ncores, pB=pB)
    OUT<-tmpSEQsw(dna_, FUN=function(x) {
          dist.xyStats(x, x.pos=x.pos_, y.pos=y.pos_, x.name=x.name,
          y.name=y.name, dist=dist, ncores=ncores)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
      end.by=end.by, wtype=wtype, global.deletion=global.deletion,
      threads=threads, pB=pB)
    return(OUT)
}
