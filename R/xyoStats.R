#' @title xyoStats
#' @name xyoStats
#' @description This function calculates \code{distIUPAC} based distances
#' comparing two populations (x: receiver and y: donor) with an outgroup
#' population (o: outgroup).
#' In the four-taxon scenario (((P1,P2),P3),O) with geneflow from P3>>P2,
#' the populations should be defined for RND, Gmin and RNDmin statistics as
#' follows [x:P2 y:P3 o:O].
#' Accordingly in the four-taxon scenario (((P1,P2),P3),O) with geneflow from
#' P2>>P3, the populations should be defined for RND, Gmin and RNDmin
#' statistics as follows [x:P3 y:P2 o:O].
#' To calculate deltaMean and deltaMin see \link[distIUPAC]{xyiStats}.
#' @param dna \code{DNAStringSet} [mandatory]
#' @param x.pos population X positions
#' [P2 population in the four-taxon scenario (((P1,P2),P3),O) with geneflow
#' from P3>>P2]  [mandatory]
#' @param y.pos population Y positions [P3 population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param o.pos population O positions [O population in the four-taxon
#' scenario (((P1,P2),P3),O) with geneflow from P3>>P2] [mandatory]
#' @param x.name population X name [default: "x"]
#' @param y.name population Y name [default: "y"]
#' @param o.name population O name [default: "o"]
#' @param chr.name chromosome name [default: "chr"]
#' @param wlen sliding windows length  [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites}  [default: "bp"]
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
#' Rosenzweig, B. K., Pease, J. B., Besansky, N. J., & Hahn, M. W. (2016).
#' Powerful methods for detecting introgressed regions from population genomic
#' data. \emph{Molecular ecology}, \bold{25(11)}, 2387-2397.
#' 
#' Paradis, E., & Schliep, K. (2018). ape 5.0: an environment for modern
#' phylogenetics and evolutionary analyses in R. \emph{Bioinformatics},
#' \bold{35(3)}, 526-528.
#'
#' Hahn, M. W., Hibbins, M. S. (2019). A Three-Sample Test for Introgression.
#' \emph{Molecular biology and evolution}, \bold{msz178}.
#' @seealso \code{\link[distIUPAC]{dist.xyoStats}},
#' \code{\link[ape]{dist.dna}}
#' @examples
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' AFG.SPRE.CAS.xyoStats<-xyoStats(MySequences, x.pos=AFG.pos,
#' y.pos=SPRE.pos, o.pos=CAS.pos, x.name="AFG", y.name = "SPRE",
#' o.name="CAS", threads=2)
#' AFG.SPRE.CAS.xyoStats
#' @export xyoStats
#' @author Kristian K Ullrich
xyoStats<-function(dna, x.pos, y.pos, o.pos, x.name="x", y.name="y",
  o.name="o", chr.name="chr", wlen=25000, wjump=25000, start.by=1, end.by=NULL,
  wtype="bp", dist="IUPAC", global.deletion=TRUE, threads=1, ncores=1,
  pB=FALSE){
    options(scipen=22)
    dna_<-dna[c(x.pos, y.pos, o.pos)]
    x.pos_<-seq(1, length(x.pos))
    y.pos_<-seq(length(x.pos_) + 1, length(x.pos_) + length(y.pos))
    o.pos_<-seq(length(x.pos_) + length(y.pos_) + 1,
      length(x.pos_) + length(y.pos_) + length(o.pos))
    #OUT<-distIUPACsw(dna_, FUN=function(x) {
    #      dist.xyoStats(x, x.pos=x.pos_, y.pos=y.pos_, o.pos=o.pos_,
    #      x.name=x.name, y.name=y.name, o.name=o.name)
    #  }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
    #  end.by=end.by, wtype=wtype, dist=dist, global.deletion=global.deletion,
    #  threads=threads, ncores=ncores, pB=pB)
    OUT<-tmpSEQsw(dna_, FUN=function(x) {
          dist.xyoStats(x, x.pos=x.pos_, y.pos=y.pos_, o.pos=o.pos_,
          x.name=x.name, y.name=y.name, o.name=o.name, dist=dist, ncores=ncores)
      }, chr.name=chr.name, wlen=wlen, wjump=wjump, start.by=start.by,
      end.by=end.by, wtype=wtype, global.deletion=global.deletion,
      threads=threads, pB=pB)
    return(OUT)
}
