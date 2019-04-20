#' @title doFasta4
#' @name doFasta4
#' @description This function uses the \code{angsd -doFasta 4} option to obtain IUPAC coded
#' fasta files from an indexed bam file and returns a \code{DNAStringSet}
#' @import Biostrings
#' @import magrittr
#' @param bam indexed bam file which needs to be indexed by 'samtools index'
#' @param bam.name individual name
#' @param angsdoption angsd option specified as here \href{angsd Filters}{http://www.popgen.dk/angsd/index.php/Filters},
#' here \href{angsd Allele Counts}{http://www.popgen.dk/angsd/index.php/Allele_Counts} and here \href{angsd doFasta}{http://www.popgen.dk/angsd/index.php/Fasta} 
#' @param region region specification ('chr' entire chr; 'chr:start' region from start to end of chr; 'chr:start-stop' region from start to stop from chr)
#' @param angsd path to angsd executable
#' @param samtools path to samtools executable
#' @param system.verbose system verbosity output 
#' @examples
#' bam.file<-system.file("extdata", "ind1.pop1.bam", package="distIUPAC")
#' #complete 
#' bam.iupac<-doFasta4(bam=bam.file, angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30 -iupacRatio 0.2", angsd="angsd", samtools="samtools")
#' #region from 5001 to 15000 from chr1 (chr1:5001-15000)
#' bam.region.iupac<-doFasta4(bam=bam.file, angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30 -iupacRatio 0.2", region="chr1:5001-15000", angsd="angsd", samtools="samtools")
#' @export doFasta4
#' @author Kristian K Ullrich
doFasta4<-function(bam, bam.name="ind1", angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30 -iupacRatio 0.2", region=NULL, angsd="angsd", samtools="samtools", system.verbose=FALSE){
  options(scipen=22)
  tmpfile<-tempfile()
  if(is.null(region)){
    if(!system.verbose){
      sprintf("%s -doFasta 4 -doCounts 1 %s -i %s -out %s", angsd, angsdoption, bam, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system(ignore.stdout=FALSE,ignore.stderr=TRUE)
      tmp.dna<-readDNAStringSet(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      names(tmp.dna)<-bam.name
      return(tmp.dna)
    }
    if(system.verbose){
      sprintf("%s -doFasta 4 -doCounts 1 %s -i %s -out %s", angsd, angsdoption, bam, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system()
      tmp.dna<-readDNAStringSet(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      names(tmp.dna)<-bam.name
      return(tmp.dna)
    }
  }
  if(!is.null(region)){
    if(!system.verbose){
      sprintf("%s -doFasta 4 -doCounts 1 %s -i %s -r %s -out %s", angsd, angsdoption, bam, region, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system(ignore.stdout=TRUE,ignore.stderr=TRUE)
      sprintf("%s faidx %s.fa.gz %s > %s", samtools, tmpfile, region, tmpfile) %>% system(ignore.stdout=FALSE,ignore.stderr=TRUE)
      tmp.dna<-readDNAStringSet(tmpfile)
      file.remove(tmpfile)
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      names(tmp.dna)<-bam.name
      return(tmp.dna)
    }
    if(system.verbose){
      sprintf("%s -doFasta 4 -doCounts 1 %s -i %s -r %s -out %s", angsd, angsdoption, bam, region, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz", samtools, tmpfile) %>% system()
      sprintf("%s faidx %s.fa.gz %s > %s", samtools, tmpfile, region, tmpfile) %>% system()
      tmp.dna<-readDNAStringSet(tmpfile)
      file.remove(tmpfile)
      file.remove(paste0(tmpfile, ".arg"))
      file.remove(paste0(tmpfile, ".fa.gz"))
      file.remove(paste0(tmpfile, ".fa.gz.fai"))
      file.remove(paste0(tmpfile, ".fa.gz.gzi"))
      names(tmp.dna)<-bam.name
      return(tmp.dna)
    }
  }
}
