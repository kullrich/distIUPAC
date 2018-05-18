#' @title faidx2seq
#' @name faidx2seq
#' @description This function returns a region of a bgzipped and samtools indexed fasta alignment file using \code{samtools}
#' @import Biostrings
#' @param fasta fasta file which needs to be indexed by 'samtools faidx'
#' @param start region start
#' @param end region end
#' @param format either \code{dna} or \code{aa}
#' @param samtools path to samtools executable
#' @examples
#' fasta.gz.file<-system.file("extdata", "example.fasta.gz", package="distIUPAC")
#' #complete alignment
#' dna<-faidx2seq(fasta.gz.file, format="dna", samtools="samtools")
#' #region from 5001 to 15000 from alignment
#' dna.region<-faidx2seq(fasta.gz.file, start=5001, end=15000, format="dna", samtools="samtools")
#' @export faidx2seq
#' @author Kristian K Ullrich
faidx2seq<-function(fasta, start=NULL, end=NULL, format="dna", samtools="samtools"){
  options(scipen=22)
  samples<-read.table(paste0(fasta,".fai"),sep="\t",header=FALSE,stringsAsFactor=FALSE)[,1]
  len<-unique(read.table(paste0(fasta,".fai"),sep="\t",header=FALSE,stringsAsFactor=FALSE)[,2])
  if(is.null(start)){start<-1}
  if(is.null(end)){end<-len}
  if(end<start){stop("end smaller than start")}
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
