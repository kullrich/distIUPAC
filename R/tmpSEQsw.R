#' @title tmpSEQsw
#' @name tmpSEQsw
#' @description This function returns tmpSEQ
#' \code{Biostrings} \code{DNAStringSet} objects
#' as sliding windows which can be used to parse to further functions
#' @import Biostrings
#' @import doMC
#' @import foreach
#' @importFrom ape dist.dna as.DNAbin
#' @importFrom stats as.dist sd setNames
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @importFrom rlist list.stack
#' @param dna \code{DNAStringSet}
#' @param FUN the function to be applied [default: NULL]
#' see e.g. \link[distIUPAC]{dist.xStats}, \link[distIUPAC]{dist.xyStats},
#' \link[distIUPAC]{dist.xyoStats}, \link[distIUPAC]{dist.xyioStats}
#' @param chr.name chromosome name [default: "chr"]
#' @param wlen sliding windows length [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' take width of \code{DNAStringSet}
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites} [default: "bp"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global or pairwise way
#' (default is to delete in a global way) [default: TRUE]
#' @param threads number of parallel threads to process windows [default: 1]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: TRUE]
#' @seealso \link[distIUPAC]{distIUPAC}, \link[ape]{dist.dna} 
#' @examples
#' data("MySequences", package="distIUPAC")
#' CAS.pos<-5:34
#' CAS.dIUPAC<-distIUPACsw(MySequences[CAS.pos])
#' ##xStats
#' CAS.xStats<-distIUPACsw(MySequences[CAS.pos], FUN=dist.xStats)
#' ##multiple threads to process windows
#' CAS.xStats<-distIUPACsw(MySequences[CAS.pos], threads=2)
#' ##multiple cores to process pairwise distance calculation
#' CAS.xStats<-distIUPACsw(MySequences[CAS.pos], ncores=2)
#' ##disbale global deletion
#' CAS.pairwiseDeletion.xStats<-distIUPACsw(MySequences[CAS.pos],
#' global.deletion=FALSE, FUN=dist.xStats, ncores=2)
#' ##using K80 distance from ape package
#' distIUPACsw(MySequences[CAS.pos], dist="K80", FUN=dist.xStats)
#' @export tmpSEQsw
#' @author Kristian K Ullrich
tmpSEQsw<-function(dna, FUN=NULL, chr.name="chr",
 wlen=25000, wjump=25000, start.by=1, end.by=NULL, wtype="bp",
 global.deletion=TRUE, threads=1, pB=TRUE){
    options(scipen=22)
    if(is.null(end.by)){end.by<-unique(width(dna))}
    if(start.by>unique(width(dna))){
        stop("start.by needs to be equal or smaller than dna length")
    }
    if(end.by>unique(width(dna))){
        stop("end.by needs to be equal or smaller than dna length")
    }
    if(wtype=="bp"){
        tmp.sw<-swgen(wlen=wlen, wjump=wjump, start.by=start.by, end.by=end.by)
    }
    if(wtype=="biSites"){
        tmp.POS<-biSites(dna, threads=threads, pB=FALSE)
        tmp.sw<-posgen(tmp.POS$biPOS, wlen=wlen, start.by=start.by,
          end.by=end.by)
    }
    if(wtype=="triSites"){
        tmp.POS<-triSites(dna, threads=threads, pB=FALSE)
        tmp.sw<-posgen(tmp.POS$triPOS, wlen=wlen, start.by=start.by,
          end.by=end.by)
    }
    j<-NULL
    if(pB){
        pb<-txtProgressBar(min=0, max=ncol(tmp.sw), initial=0, style=3)
    }
    registerDoMC(threads)
    OUT<-foreach(j=seq(from=1, to=ncol(tmp.sw)), .combine=rbind) %dopar% {
        CHRNAME<-chr.name
        START<-tmp.sw[1, j][[1]]
        END<-tmp.sw[2, j][[1]]
        distIUPAC<-NA
        OUT<-list(CHRNAME, START, END, distIUPAC)
        names(OUT)<-c("CHRNAME", "START", "END", "tmpSEQ")
        tmp.seq<-subseq(dna, START, END)
        if(global.deletion){
            tmp.seq<-globalDeletion(tmp.seq, pB=FALSE)
        }
        OUT$tmpSEQ<-tmp.seq
        if(pB){
            setTxtProgressBar(pb, j)
        }
        OUT
    }
    if(pB){
        setTxtProgressBar(pb, ncol(tmp.sw))
        close(pb)
    }
    if(is.null(FUN)){
        return(OUT)
    }
    if(is.function(FUN)){
        if(!is.matrix(OUT)){
            return(
              c(OUT[1:3], rlist::list.stack(lapply(OUT[4], FUN)))
            )
        }
        return(
          cbind(OUT[, 1:3], rlist::list.stack(lapply(OUT[, 4], FUN)))
        )
    }
}
