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
- [RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)
- [ape](https://cran.r-project.org/web/packages/ape/index.html)
- [foreach](https://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](https://cran.r-project.org/web/packages/doMC/index.html)

```
install.packages("devtools")
install.packages("Rcpp")
install.packages("RcppParallel")
install.packages("ape")
install.packages("foreach")
install.packages("doMC")
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
- [02. distStats: population distanceStats per sliding-window: iupac fasta + pop info >>> dStats](https://github.com/kullrich/distIUPAC/tree/master/vignettes/dStats.Rmd)
- [03. getTrees: genrate twisst input: iupac fasta >>> twisst trees](https://github.com/kullrich/distIUPAC/tree/master/vignettes/twisstTrees.Rmd)
- [04. getTreemix: genrate treemix input: iupac fasta >>> treemix](https://github.com/kullrich/distIUPAC/tree/master/vignettes/treemix.Rmd)
- [05. distFeatures: using GFF3/GTF featues: iupac fasta + GTF >>> feature specific distances](https://github.com/kullrich/distIUPAC/tree/master/vignettes/GTFdistances.Rmd)

### Quick-guide

```
library(distIUPAC)

#browse vignettes
browseVignettes("distIUPAC")

#load IUPAC encoded nucleotide sequences with Biostrings
#change path to your input file
input.fasta <- paste0(find.package("distIUPAC"),"/data/seqIUPAC.fasta")
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
xStats(MySequences[CAS.pos])

#get mean IUPAC distances using a pre-defined distance matrix on sliding windows using multiple threads (xStats)
xStats(MySequences[CAS.pos], threads = 4)

#get bi-allelic sites for population
CAS.biSites<-biSites(MySequences, x.pos = CAS.pos)
as.matrix(MySequences[CAS.pos])[,head(CAS.biSites)]
```

### License

MIT (see LICENSE)

### Bug reports

Please report any errors or requests regarding `distIUPAC` to Kristian Ullrich (ullrich@evolbio.mpg.de)
