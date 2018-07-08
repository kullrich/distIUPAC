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
#' @examples
#' data("MySequences", package = "distIUPAC")
#' CAS.pos<-5:34
#' AFG.pos<-82:87
#' SPRE.pos<-106:113
#' pop.list<-list(CAS.pos, AFG.pos, SPRE.pos)
#' names(pop.list)<-c("CAS", "AFG", "SPRE")
#' CAS.AFG.SPRE.treemix<-iupac2treemix(MySequences, list.pos = pop.list,
#' wlen = 500, threads = 2)
#' head(CAS.AFG.SPRE.treemix)
#' @export iupac2treemix
#' @author Kristian K Ullrich
iupac2treemix<-function(dna, list.pos, wlen=25000, start.by=NULL, end.by=NULL, threads=1){
  IUPAC_CODE_MAP_LIST<-list(c("A","A"),c("C","C"),c("G","G"),c("T","T"),c("A","C"),c("A","G"),c("A","T"),c("C","G"),c("C","T"),c("G","T"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"),c("N","N"))
  names(IUPAC_CODE_MAP_LIST)<-c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B","N","-","+",".")
  options(scipen=22)
  CONVERT<-function(x, list.pos){
    j<-rev(x)[[1]]
    biPOS<-rev(x)[[2]]
    biA<-NA
    biB<-NA
    biFRQ<-rep(list(NA),length(list.pos))
    names(biFRQ)<-names(list.pos)
    OUT<-append(list(biPOS,biA,biB),biFRQ)
    names(OUT)<-c("biPOS","biA","biB",names(biFRQ))
    biAB<-names(table(unlist(lapply(x[unlist(list.pos)],function(y) IUPAC_CODE_MAP_LIST[y]))))
    biAB<-biAB[biAB%in%c("A","C","G","T")]
    OUT$biA<-biAB[1]
    OUT$biB<-biAB[2]
    pop.biAB<-lapply(list.pos,function(y) table(unlist(lapply(x[y],function(z) IUPAC_CODE_MAP_LIST[z]))))
    pop.biA<-lapply(pop.biAB,function(y) max(ifelse(names(y)==OUT$biA,y,0)))
    pop.biB<-lapply(pop.biAB,function(y) max(ifelse(names(y)==OUT$biB,y,0)))
    pop.out<-paste(pop.biA,pop.biB,sep=",")
    names(pop.out)<-names(list.pos)
    OUT[names(pop.out)]<-pop.out
    OUT
  }
  if(is.null(names(list.pos))){names(list.pos)<-seq(1,length(list.pos))}
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  dna.biSites<-biSites(dna[unlist(list.pos)], wlen=wlen, start.by=start.by, end.by=end.by, threads=threads, pB=FALSE)
  dna_<-rbind(as.matrix(subString(dna,dna.biSites,dna.biSites)),dna.biSites,1:length(dna.biSites))
  cl<-makeCluster(threads)
  clusterExport(cl, list("dna_","list.pos","CONVERT","IUPAC_CODE_MAP_LIST"), envir=environment())
  dna.treemix<-parApply(dna_, 2, function(x) CONVERT(x=x, list.pos=list.pos), cl=cl)
  stopCluster(cl)
  OUT<-foreach(j=1:length(dna.treemix), .combine=rbind) %do% {
    unlist(dna.treemix[j])
  }
  return(OUT)
}
