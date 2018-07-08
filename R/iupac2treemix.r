#' @title iupac2treemix
#' @name iupac2treemix
#' @description This function returns biallelic sites as treemix-like output.
#' @import Biostrings
#' @import ape
#' @import parallel
#' @import foreach
#' @param dna \code{DNAStringSet}
#' @param list.pos population positions list
#' @param wlen sliding window length for biSites retrieval
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
#' CAS.AFG.SPRE.treemix<-iupac2treemix(MySequences, list.pos = pop.list, threads = 2)
#' head(CAS.AFG.SPRE.treemix)
#' @export iupac2treemix
#' @author Kristian K Ullrich
iupac2treemix<-function(dna, list.pos, wlen=25000, start.by=NULL, end.by=NULL, threads=1, pB=TRUE){
  IUPAC_CODE_MAP_LIST<-list(c("A"),c("C"),c("G"),c("T"),c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"),c("A","C","G"),c("A","C","T"),c("A","G","T"),c("C","G","T"),c(),c(),c(),c())
  names(IUPAC_CODE_MAP_LIST)<-c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B","N","-","+",".")
  IUPAC_CODE_MAP_LIST_<-list(c("A","A"),c("C","C"),c("G","G"),c("T","T"),c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"))
  names(IUPAC_CODE_MAP_LIST_)<-c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B","N","-","+",".")
  options(scipen=22)
  if(is.null(names(list.pos))){names(list.pos)<-seq(1,length(list.pos))}
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
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
    START<-NA
    END<-NA
    biPOS<-NA
    OUT<-list(START,END,biPOS)
    names(OUT)<-c("START","END","biPOS")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(unoque(width(tmp.seq))==1){
      tmp.seq.cM<-t(as.matrix(apply(consensusMatrix(tmp.seq),1,function(x) ifelse(x>0,1,0))))
    }
    if(unoque(width(tmp.seq))!=1){
      tmp.seq.cM<-apply(consensusMatrix(tmp.seq),1,function(x) ifelse(x>0,1,0))
    }
    tmp.biPOS<-OUT$START-1+which(apply(tmp.seq.cM,1,function(x) length(unique(unlist(unique(IUPAC_CODE_MAP_LIST[names(x[x==1])])))))==2)
    if(length(tmp.biPOS)==0){return(NULL)}
    tmp.biPOS_<-which(apply(tmp.seq.cM,1,function(x) length(unique(unlist(unique(IUPAC_CODE_MAP_LIST[names(x[x==1])])))))==2)
    tmp.seq.biPOS<-subString(tmp.seq,tmp.biPOS_,tmp.biPOS_)
    k<-NULL
    TREEMIX<-foreach(k=1:length(tmp.biPOS_), .combine=rbind) %do% {
      tmp.bi<-subseq(tmp.seq.biPOS,k,k)
      biPOS<-tmp.biPOS[k]
      biA<-NA
      biB<-NA
      biFRQ<-rep(list(NA),length(pos_))
      names(biFRQ)<-names(pos_)
      TREEMIX<-append(list(biPOS,biA,biB),biFRQ)
      names(TREEMIX)<-c("biPOS","biA","biB",names(biFRQ))
      biAB<-names(table(unlist(lapply(as.character(tmp.bi),function(y) IUPAC_CODE_MAP_LIST_[y]))))
      biAB<-biAB[biAB%in%c("A","C","G","T")]
      TREEMIX$biA<-biAB[1]
      TREEMIX$biB<-biAB[2]
      pop.biAB<-lapply(pos_,function(y) table(unlist(lapply(as.character(tmp.bi)[y],function(z) IUPAC_CODE_MAP_LIST_[z]))))
      pop.biA<-lapply(pop.biAB,function(y) max(ifelse(names(y)==TREEMIX$biA,y,0)))
      pop.biB<-lapply(pop.biAB,function(y) max(ifelse(names(y)==TREEMIX$biB,y,0)))
      pop.out<-paste(pop.biA,pop.biB,sep=",")
      names(pop.out)<-names(pos_)
      TREEMIX[names(pop.out)]<-pop.out
      TREEMIX
    }
    if(pB){
      setTxtProgressBar(pb,j)      
    }
    TREEMIX
  }
  if(pB){
    setTxtProgressBar(pb,dim(tmp.sw)[2])
    close(pb)    
  }
  return(OUT)
}
