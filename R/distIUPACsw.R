#' @title distIUPACsw
#' @name distIUPACsw
#' @description This function calculates \code{IUPAC} distances on \code{Biostrings} \code{DNAStringSet} objects
#' on sliding windows which can be used to parse to further functions
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @import magrittr
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @importFrom rlist list.stack
#' @param dna \code{DNAStringSet}
#' @param FUN the function to be applied
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param start.by optional start position
#' @param end.by optional end position
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites with missing data in a global or pairwise way (default is to delete in a global way)
#' @param threads number of parallel threads
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' @export distIUPACsw
#' @author Kristian K Ullrich
distIUPACsw<-function(dna, FUN=NULL, chr.name="chr",
 wlen=25000, wjump=25000, start.by=NULL, end.by=NULL, wtype="bp",
 dist="IUPAC", global.deletion=TRUE, threads=1, pB=TRUE){
  options(scipen=22)
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=start.by,end.by=end.by)
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$biPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS$triPOS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  j<-NULL
  if(pB){
    pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)
  }
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
    CHRNAME<-chr.name
    START<-tmp.sw[1,j][[1]]
    END<-tmp.sw[2,j][[1]]
    distIUPAC<-NA
    OUT<-list(CHRNAME,START,END,distIUPAC)
    names(OUT)<-c("CHRNAME","START","END","distIUPAC")
    tmp.seq<-subseq(dna,START,END)
    if(global.deletion){
      tmp.seq<-globalDeletion(tmp.seq,pB=FALSE)
    }
    if(dist=="IUPAC"){
      OUT$distIUPAC<-distIUPAC(as.character(tmp.seq))
    }
    if(dist!="IUPAC"){
      OUT$distIUPAC<-setNames( list(dist.dna(as.DNAbin(DNAMultipleAlignment(tmp.seq)),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE),
       pairwiseDeletion(as.character(tmp.seq))$sitesUsed), c("distIUPAC","sitesUsed"))
    }
    if(pB){
      setTxtProgressBar(pb,j)
    }
    OUT
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)
  }
  if(is.null(FUN)){
    return(OUT)
  }
  if(is.function(FUN)){
  	if(!is.matrix(OUT)){
  	  return(c(OUT[1:3],OUT[4] %>% lapply(., FUN) %>% rlist::list.stack(.)))	
  	}
    return(cbind(OUT[,1:3],OUT[,4] %>% lapply(., FUN) %>% rlist::list.stack(.)))
  }
}
