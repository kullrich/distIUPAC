#' @title faidx2seq
#' @name faidx2seq
#' @description This function returns a region of a bgzipped and samtools indexed fasta alignment file using \code{samtools}
#' @import Biostrings
#' @param fasta fasta file which needs to be indexed by 'samtools faidx'
#' @param start region start
#' @param end region end
#' @param format either dna or aa
#' @examples
#' @export faidx2seq
#' @author Kristian K Ullrich
faidx2seq<-function(fasta, start, end, format="dna", samtools="samtools"){
  options(scipen=22)
  if(end>start){stop("end smaller than start")}
  samples<-read.table(paste0(fasta,".fai"),sep="\t",header=FALSE,stringsAsFactor=FALSE)[,1]
  tmpfile<-tempfile()
  for(i in samples){
    system(paste0(samtools," faidx ",fasta," ",i,":",start,"-",end," >> ",tmpfile))
  }
  if(format=="dna"){
    out<-readDNAStringSet(tmpfile)
  }
  if(format=="aa"){
    out<-readAAStringSet(tmpfile)
  }
  system(paste0("rm ",tmpfile))
  return(out)
}
