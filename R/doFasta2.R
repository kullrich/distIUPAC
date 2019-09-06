#' @title doFasta2
#' @name doFasta2
#' @description This function uses the \code{angsd -doFasta 2} option
#' to obtain a consensus fasta file from a bam file list of indexed bam files
#' and returns a \code{DNAStringSet}.
#' @import Biostrings
#' @param bamlist list of indexed bam files which
#' needs to be indexed by 'samtools index'
#' @param bamlist.name population name [default: "pop1"]
#' @param angsdoption angsd option specified as follows
#' \href{angsd Filters}{http://www.popgen.dk/angsd/index.php/Filters},
#' \href{Allele Counts}{http://www.popgen.dk/angsd/index.php/Allele_Counts}
#' and \href{angsd doFasta}{http://www.popgen.dk/angsd/index.php/Fasta}
#' [default: "-maxDepth 99999 -minQ 20 -minMapQ 30"]
#' @param region region specification ('chr' entire chr; 'chr:start'
#' region from start to end of chr; 'chr:start-stop'
#' region from start to stop from chr) [default: NULL]
#' @param angsd path to angsd executable [default: "angsd"]
#' @param samtools path to samtools executable [default: "samtools"]
#' @param system.verbose system verbosity output [default: FALSE]
#' @examples
#' bam.url<-"http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/bam/"
#' bam.ind1<-"SS6004467-dedup.rg.bam"
#' bam.ind2<-"SS6004468-dedup.rg.bam"
#' bam.files<-c(paste0(bam.url, bam.ind1), paste0(bam.url, bam.ind2))
#' ##region from 5001 to 15000 from chr1 (1:5001-15000)
#' #bam.region.consensus<-doFasta2(bamlist=bam.files,
#' #  angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30",
#' #  region="1:5001-15000", angsd="angsd", samtools="samtools")
#' @export doFasta2
#' @author Kristian K Ullrich
doFasta2<-function(bamlist, bamlist.name="pop1",
  angsdoption="-maxDepth 99999 -minQ 20 -minMapQ 30", region=NULL,
  angsd="angsd", samtools="samtools", system.verbose=FALSE){
    options(scipen=22)
    if(!file.exists(angsd)){
        stop("angsd not in path")
    }
    if(!file.exists(samtools)){
        stop("samtools not in path")
    }
    tmpfile<-tempfile()
    bamlistfile<-tempfile()
    cat(bamlist, sep="\n", file=bamlistfile)
    if(is.null(region)){
        angsd.args<-sprintf("-doFasta 2 -doCounts 1 %s -b %s -out %s",
          angsdoption, bamlistfile, tmpfile)
        samtools.args<-sprintf("faidx %s.fa.gz", tmpfile)
    }
    if(!is.null(region)){
        angsd.args<-sprintf("-doFasta 2 -doCounts 1 %s -b %s -r %s -out %s",
          angsdoption, bamlistfile, region, tmpfile)
        samtools.args1<-sprintf("faidx %s.fa.gz", tmpfile)
        samtools.args2<-sprintf("faidx %s.fa.gz %s > %s", tmpfile, region,
          tmpfile)
    }
    if(!system.verbose && is.null(region)){
        system2(command=angsd, args=angsd.args, stdout=NULL, stderr=NULL)
        system2(command=samtools, args=samtools.args, stderr=NULL)
    }
    if(system.verbose && is.null(region)){
        system2(command=angsd, args=angsd.args)
        system2(command=samtools, args=samtools.args)
    }
    if(!system.verbose && !is.null(region)){
        system2(command=angsd, args=angsd.args, stdout=NULL, stderr=NULL)
        system2(command=samtools, args=samtools.args1, stdout=NULL,
          stderr=NULL)
        system2(command=samtools, args=samtools.args2, stderr=NULL)
    }
    if(system.verbose && !is.null(region)){
        system2(command=angsd, args=angsd.args)
        system2(command=samtools, args=samtools.args1)
        system2(command=samtools, args=samtools.args2)
    }
    if(is.null(region)){
        tmp.dna<-readDNAStringSet(paste0(tmpfile, ".fa.gz"))
    }
    if(!is.null(region)){
        tmp.dna<-readDNAStringSet(tmpfile)
    }
    file.remove(paste0(tmpfile, ".arg"))
    file.remove(paste0(tmpfile, ".fa.gz"))
    file.remove(paste0(tmpfile, ".fa.gz.fai"))
    file.remove(paste0(tmpfile, ".fa.gz.gzi"))
    file.remove(tmpfile)
    file.remove(bamlistfile)
    names(tmp.dna)<-bamlist.name
    return(tmp.dna)
}
