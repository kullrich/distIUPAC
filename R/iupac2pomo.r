#' @title iupac2pomo
#' @name iupac2pomo
#' @description This function returns sites as IQ-TREE PoMo-like output.
#' @import Biostrings
#' @import ape
#' @import parallel
#' @import foreach
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param list.pos population positions list
#' @param chr set chromosome name
#' @param wlen sliding window length for PoMo sites retrieval
#' @param start.by optional start position
#' @param end.by optional end position
#' @param threads number of parallel threads
#' @param pB specifies if progress should be shown as a progress bar
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' pop.list<-list(CAS.pos, AFG.pos, SPRE.pos)
#' names(pop.list)<-c("CAS", "AFG", "SPRE")
#' CAS.AFG.SPRE.pomo<-iupac2pomo(MySequences, list.pos = pop.list, threads = 2)
#' head(CAS.AFG.SPRE.pomo)
#' @export iupac2pomo
#' @author Kristian K Ullrich
iupac2pomo<-function(dna, list.pos, chr="1", wlen=25000, start.by=NULL, end.by=NULL, threads=1, pB=TRUE){
  diploid2haploid<-function(x){
    x["A"]<-x["A"]*2
    x["C"]<-x["C"]*2
    x["G"]<-x["G"]*2
    x["T"]<-x["T"]*2
    x["N"]<-x["N"]+x["-"]
    x["-"]<-0
    x["N"]<-x["N"]+x["+"]
    x["+"]<-0
    x["N"]<-x["N"]+x["."]
    x["."]<-0
    x["N"]<-x["N"]*2
    x["A"]<-x["A"]+x["M"]
    x["C"]<-x["C"]+x["M"]
    x["M"]<-0;
    x["A"]<-x["A"]+x["R"]
    x["G"]<-x["G"]+x["R"]
    x["R"]<-0
    x["A"]<-x["A"]+x["W"]
    x["T"]<-x["T"]+x["W"]
    x["W"]<-0
    x["C"]<-x["C"]+x["S"]
    x["G"]<-x["G"]+x["S"]
    x["S"]<-0;
    x["C"]<-x["C"]+x["Y"]
    x["T"]<-x["T"]+x["Y"]
    x["Y"]<-0
    x["G"]<-x["G"]+x["K"]
    x["T"]<-x["T"]+x["K"]
    x["K"]<-0
    x["N"]<-x["N"]+(x["V"]*2)
    x["V"]<-0
    x["N"]<-x["N"]+(x["H"]*2)
    x["H"]<-0
    x["N"]<-x["N"]+(x["D"]*2)
    x["D"]<-0
    x["N"]<-x["N"]+(x["B"]*2)
    x["B"]<-0
    return(x)
  }
  options(scipen=22)
  if(is.null(names(list.pos))){names(list.pos)<-seq(1,length(list.pos))}
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  dna_<-dna[unlist(list.pos)]
  list.pos_<-unlist(lapply(list.pos,length))
  cur_<-1
  pos_<-list()
  for(i in 1:length(list.pos_)){
    if(i==1){
      pos_<-append(pos_,list(seq(from=cur_,to=list.pos_[i])))
      cur_<-list.pos_[i]+1
    }
    if(i!=1){
      pos_<-append(pos_,list(seq(from=cur_,to=cur_+list.pos_[i]-1)))
      cur_<-cur_+list.pos_[i]
    }
  }
  names(pos_)<-names(list.pos)
  tmp.sw<-swgen(wlen=wlen,wjump=wlen,start.by=start.by,end.by=end.by)
  if(pB){
    pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)  
  }
  j<-NULL
  registerDoMC(threads)
  OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
    CHROM<-chr
    START<-NA
    END<-NA
    OUT<-list(CHROM,START,END)
    names(OUT)<-c("CHROM","START","END")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    pop.cM<-lapply(pos_,function(y) apply(consensusMatrix(tmp.seq[y]),2,diploid2haploid)[1:4,])
    k<-NULL
    POMO<-foreach(k=1:length(pop.cM), .combine=cbind) %do% {
      POMO<-apply(pop.cM[[k]],2,function(y) paste(y,collapse=","))
      POMO
    }
    colnames(POMO)<-names(pos_)
    OUT<-cbind(OUT$CHROM,seq(from=OUT$START,to=OUT$END),POMO)
    colnames(OUT)[1]<-"CHROM"
    colnames(OUT)[2]<-"POS"
    if(pB){
      setTxtProgressBar(pb,j)      
    }
    OUT
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)    
  }
  return(OUT)
}
