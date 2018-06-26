distIUPAC
=========
### 

### Installation

## R specific installation prerequisites

First install [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) package from [bioconductor](https://bioconductor.org/).

```
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
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
#install_github("kullrich/distIUPAC", build_vignettes = TRUE, dependencies = FALSE)
install_git("https://gwdg.gitlab.de/evolgen/distIUPAC.git", build_vignettes = TRUE, dependencies = FALSE)
```

### Vignettes

These vignettes introduce `distIUPAC`

- [00. Intro](https://github.com/kullrich/distIUPAC/tree/master/vignettes/Intro.Rmd)
- [01. bam2IUPAC: fastq + reference >>> reference mapped bam + angsd >>> iupac fasta](https://github.com/kullrich/distIUPAC/tree/master/vignettes/bam2IUPAC.Rmd)
- [02. distStats: population distanceStats per sliding-window: iupac fasta + pop info >>> dStats ](https://github.com/kullrich/distIUPAC/tree/master/vignettes/dStats.Rmd)
- [03. getTrees: genrate twisst input: iupac fasta >>> twisst trees](https://github.com/kullrich/distIUPAC/tree/master/vignettes/twisstTrees.Rmd)
- [04. distFeatures: using GFF3/GTF featues: iupac fasta + GTF >>> feature specific distances ](https://github.com/kullrich/distIUPAC/tree/master/vignettes/GTFdistances.Rmd)

### Quick-guide

```
library(distIUPAC)

#browse vignettes
browseVignettes("distIUPAC")

#load IUPAC encoded nucleotide sequences with Biostrings
#change path to your input file
input.fasta <- paste0(find.package("distIUPAC"),"/data/seqIUPAC.fasta")
MySequences <- readBStringSet(input.fasta)
MySequences

#get IUPAC distances using a pre-defined distance matrix
distIUPAC(as.character(MySequences))

#get IUPAC distance using your own distance matrix
MyScoreMatrix <- scoreMatrix()
MyScoreMatrix["C","Y"] <- 1.0
distIUPACmatrix(as.character(MySequences), MyScoreMatrix)
```

### License

MIT (see LICENSE)

### Bug reports

Please report any errors or requests regarding `distIUPAC` to Kristian Ullrich (ullrich@evolbio.mpg.de)
