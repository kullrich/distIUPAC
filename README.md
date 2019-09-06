distIUPAC
=========
### 

### Installation

## R specific installation prerequisites

First install [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package from [bioconductor](https://bioconductor.org/).

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings", version = "3.8")
```

Second install additional packages from cran [CRAN](https://cran.r-project.org/web/packages/index.html).
- [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [RcppThread](https://cran.r-project.org/web/packages/RcppThread/index.html)
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](https://cran.r-project.org/web/packages/doMC/index.html)

```
install.packages("devtools")
install.packages("Rcpp")
install.packages("RcppThread")
install.packages("ape")
install.packages("foreach")
install.packages("doMC")
install.packages("magrittr")
install.packages("rlist")
```

Third install `distIUPAC` package from [github](https://github.com/kullrich) or [gwdg gitlab](https://gwdg.gitlab.de) using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```
library(devtools)
install_github("kullrich/distIUPAC", build_vignettes = TRUE, dependencies = FALSE)
#install_git("https://gwdg.gitlab.de/evolgen/distIUPAC.git", build_vignettes = TRUE, dependencies = FALSE)
```

### Vignettes

These vignettes introduce `distIUPAC`

- [00. Intro](https://github.com/kullrich/distIUPAC/tree/master/vignettes/Intro.Rmd)
- [01. bam2IUPAC: fastq + reference >>> reference mapped bam + angsd >>> iupac fasta](https://github.com/kullrich/distIUPAC/tree/master/vignettes/bam2IUPAC.Rmd)
- [02. distStats+ABBA-BABA: population distanceStats per sliding-window: iupac fasta + pop info >>> dStats](https://github.com/kullrich/distIUPAC/tree/master/vignettes/dStats.Rmd)
- [03. getTrees: genrate twisst input: iupac fasta >>> twisst trees](https://github.com/kullrich/distIUPAC/tree/master/vignettes/twisstTrees.Rmd)
- [04. getTreemix: genrate treemix input: iupac fasta >>> treemix](https://github.com/kullrich/distIUPAC/tree/master/vignettes/treemix.Rmd)
- [05. distFeatures: using GFF3/GTF features: iupac fasta + GTF >>> feature specific distances](https://github.com/kullrich/distIUPAC/tree/master/vignettes/GTFdistances.Rmd)

### Quick-guide

```
library(distIUPAC)

#browse vignettes
browseVignettes("distIUPAC")

#load IUPAC encoded nucleotide sequences with Biostrings
#change path to your input file
#lada fasta file
input.fasta <- readDNAStringSet(paste0(find.package("distIUPAC"),"/data/seqIUPAC.fasta.gz"))


#use example sequences
data("MySequences", package = "distIUPAC")
MySequences

#consider only a subset of all sequences
CAS.pos <- 5:34

#get IUPAC distances using a pre-defined distance matrix
CAS.distIUPAC <- distIUPAC(as.character(MySequences[CAS.pos]))

#get pairwise IUPAC distances as distance matrix
as.dist(CAS.distIUPAC$distIUPAC)

#get pairwise used sites as distance matrix
as.dist(CAS.distIUPAC$sitesUsed)

#plot bionj tree based on IUPAC distances
plot(bionj(as.dist(CAS.distIUPAC$distIUPAC)))

#get IUPAC distance using your own distance matrix
MyScoreMatrix <- scoreMatrix()
MyScoreMatrix["C","Y"] <- 1.0
distIUPACmatrix(as.character(MySequences[CAS.pos]), MyScoreMatrix)

#get mean IUPAC distances using a pre-defined distance matrix on sliding windows (xStats)
xStats(MySequences, x.pos=CAS.pos)

#get mean IUPAC distances using a pre-defined distance matrix on sliding windows using multiple threads (xStats)
xStats(MySequences[CAS.pos], threads = 4)

#get dXY calculations for three populations ((P1:i,P2:x),P3:y) on sliding windows and the (xyiStats)
AFG.pos<-82:87
SPRE.pos<-106:113
xyiStats(MySequences, x.pos=AFG.pos, y.pos=SPRE.pos, i.pos=CAS.pos, threads = 4)

#calculate ABBA-BABA statistics for four-taxon scenario (((P1:i,P2:x),P3:y),P4:o)
APO.pos<-1
xyiStats(MySequences, x.pos=AFG.pos, y.pos=SPRE.pos, i.pos=CAS.pos, o.pos=APO.pos threads = 4, do.ABBA = TRUE)

#get bi-allelic sites for population
CAS.biSites<-biSites(MySequences, x.pos = CAS.pos)
as.matrix(MySequences[CAS.pos])[,head(CAS.biSites)]
```

### Todo
- add mask file for sliding window calculations

### License

MIT (see LICENSE)

### Bug reports

Please report any errors or requests regarding `distIUPAC` to Kristian Ullrich (ullrich@evolbio.mpg.de)
