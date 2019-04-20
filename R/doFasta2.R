#' @title doFasta2
#' @name doFasta2
#' @description This function uses the \code{angsd -doFasta 2} option to obtain a consensus
#' fasta file from an bam file list of indexed bam files and returns a \code{DNAStringSet}
#' @import Biostrings
#' @import magrittr
#' @param bamlist list of indexed bam files which needs to be indexed by 'samtools index'
#' @param bamlist.name population name
#' @param angsdoption angsd option specified as here \href{angsd Filters}{http://www.popgen.dk/angsd/index.php/Filters},
#' here \href{angsd Allele Counts}{http://www.popgen.dk/angsd/index.php/Allele_Counts} and here \href{angsd doFasta}{http://www.popgen.dk/angsd/index.php/Fasta} 
#' @param region region specification ('chr' entire chr; 'chr:start' region from start to end of chr; 'chr:start-stop' region from start to stop from chr)
#' @param angsd path to angsd executable
#' @param samtools path to samtools executable
#' @param system.verbose system verbosity output 
#' @examples
#' bam.list<-c(system.file("extdata", "ind1.pop1.bam", package="distIUPAC"),system.file("extdata", "ind2.pop1.bam", package="distIUPAC"))
#' #complete 
#' bam.consensus<-doFasta2(bamlist=bam.list, angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30", angsd="angsd", samtools="samtools")
#' #region from 5001 to 15000 from chr1 (chr1:5001-15000)
#' bam.region.consensus<-doFasta2(bamlist=bam.list, angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30", region="chr1:5001-15000", angsd="angsd", samtools="samtools")
#' @export doFasta2
#' @author Kristian K Ullrich
doFasta2<-function(bamlist, bamlist.name="pop1", angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30", region=NULL, angsd="angsd", samtools="samtools", system.verbose=FALSE){
  options(scipen=22)
  tmpfile<-tempfile()
  bamlistfile<-tempfile()
  cat(bamlist, sep="\n", file=bamlistfile)
  if(is.null(region)){
    if(!system.verbose){
      sprintf("%s -doFasta 2 -doCounts 1 %s -b %s -out %s", angsd, angsdoption, bamlistfile, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system(ignore.stdout=FALSE,ignore.stderr=TRUE)
      tmp.dna<-readDNAStringSet(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      file.remove(bamlistfile)
      names(tmp.dna)<-bamlist.name
      return(tmp.dna)
    }
    if(system.verbose){
      sprintf("%s -doFasta 2 -doCounts 1 %s -b %s -out %s", angsd, angsdoption, bamlistfile, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system()
      tmp.dna<-readDNAStringSet(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      file.remove(bamlistfile)
      names(tmp.dna)<-bamlist.name
      return(tmp.dna)
    }
  }
  if(!is.null(region)){
    if(!system.verbose){
      sprintf("%s -doFasta 2 -doCounts 1 %s -b %s -r %s -out %s", angsd, angsdoption, bamlistfile, region, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz %s > %s", samtools, tmpfile, region, tmpfile) %>% system(ignore.stdout=FALSE,ignore.stderr=TRUE)
      tmp.dna<-readDNAStringSet(tmpfile)
      file.remove(tmpfile)
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      file.remove(bamlistfile)
      names(tmp.dna)<-bamlist.name
      return(tmp.dna)
    }
    if(system.verbose){
      sprintf("%s -doFasta 2 -doCounts 1 %s -b %s -r %s -out %s", angsd, angsdoption, bamlistfile, region, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz %s > %s", samtools, tmpfile, region, tmpfile) %>% system()
      tmp.dna<-readDNAStringSet(tmpfile)
      file.remove(tmpfile)
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      file.remove(bamlistfile)
      names(tmp.dna)<-bamlist.name
      return(tmp.dna)
    }
  }
}
