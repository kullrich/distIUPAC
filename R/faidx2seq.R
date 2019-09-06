#' @title faidx2seq
#' @name faidx2seq
#' @description This function returns a region of a bgzipped and samtools
#' indexed fasta alignment file using \code{samtools}
#' @import Biostrings
#' @param fasta fasta file which needs to be indexed by 'samtools faidx'
#' @param start region start [default: 1]
#' @param end region end [default: NULL]
#' @param format either \code{dna} or \code{aa} [default: "dna"]
#' @param samtools path to samtools executable [default: "samtools"]
#' @examples
#' fasta.gz.file<-system.file("extdata", "example.fasta.gz", 
#'   package="distIUPAC")
#' ##complete alignment
#' dna<-faidx2seq(fasta.gz.file, format="dna", samtools="samtools")
#' ##region from 5001 to 15000 from alignment
#' dna.region<-faidx2seq(fasta.gz.file, start = 5001, end = 15000,
#'   format = "dna", samtools = "samtools")
#' @export faidx2seq
#' @author Kristian K Ullrich
faidx2seq<-function(fasta, start=1, end=NULL, format="dna",
  samtools="samtools"){
    options(scipen=22)
    samples<-read.table(paste0(fasta, ".fai"), sep="\t", header=FALSE,
      stringsAsFactors=FALSE)[, 1]
    len<-unique(read.table(paste0(fasta, ".fai"), sep="\t", header=FALSE,
      stringsAsFactors=FALSE)[, 2])
    if(is.null(end)){end<-len}
    if(end<start){stop("end smaller than start")}
    tmpfile<-tempfile()
    for(s in samples){
        samtools.args<-sprintf("faidx %s %s:%s-%s >> %s", fasta, s, start, end,
          tmpfile)
        system2(command=samtools, args=samtools.args)
    }
    if(format=="dna"){
        out<-readDNAStringSet(tmpfile)
    }
    if(format=="aa"){
        out<-readAAStringSet(tmpfile)
    }
    file.remove(tmpfile)
    return(out)
}
