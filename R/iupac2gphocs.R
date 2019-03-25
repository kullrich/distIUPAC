#' @title iupac2gphocs
#' @name iupac2gphocs
#' @description This function returns loci blocks as G-PhoCS-like output.
#' @import Biostrings
#' @import ape
#' @import parallel
#' @import foreach
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param chr set chromosome name
#' @param wlen sliding window length for loci retrieval
#' @param start.by optional start position
#' @param end.by optional end position
#' @param threads number of parallel threads
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' gphocs<-iupac2gphocs(MySequences, threads = 2)
#' str(gphocs)
#' #sink("seqs-sample.txt")
#' #cat(gphocs,sep="\n")
#' #sink(NULL)
#' @export iupac2gphocs
#' @author Kristian K Ullrich
iupac2gphocs<-function(dna, chr="1", wlen=25000, start.by=NULL, end.by=NULL, threads=1, pB=TRUE){
  options(scipen=22)
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  dna_<-dna
  tmp.sw<-swgen(wlen=wlen,wjump=wlen,start.by=start.by,end.by=end.by)
  if(pB){
    pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)
  }
  j<-NULL
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=c) %dopar% {
    CHROM<-chr
    START<-NA
    END<-NA
    LOCUSID<-NA
    OUT<-list(CHROM,START,END,LOCUSID)
    names(OUT)<-c("CHROM","START","END","LOCUSID")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    OUT$LOCUSID<-paste0(OUT$CHROM,"_",OUT$START,"_",OUT$END," ",length(dna_)," ",OUT$END-OUT$START+1)
    tmp.seq.char<-as.character(subseq(dna_,OUT$START,OUT$END))
    LOCUS<-paste0(names(tmp.seq.char)," ",tmp.seq.char)
    OUT<-c(OUT$LOCUSID,LOCUS)
    if(pB){
      setTxtProgressBar(pb,j)      
    }
    OUT
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)    
  }
  return(c(dim(tmp.sw)[2],OUT))
}
