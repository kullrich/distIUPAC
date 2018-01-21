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

Second install [ape](http://ape-package.ird.fr/) package.

```
install.packages("ape")
```

Third install `distIUPAC` package from [github](https://github.com/kullrich) using the [devtools](https://cran.r-project.org/web/packages/devtools/index.html) package.

```
#install.packages("devtools")
library(devtools)
install_github("kullrich/distIUPAC", build_vignettes = TRUE, dependencies = TRUE)
```

During the installation process following R packages will be installed as dependencies:

* [Rcpp](Rcpp)
* [RcppParallel]()

### Vignettes

These vignettes introduce `distIUPAC`

- [distIUPAC Intro](https://github.com/kullrich/distIUPAC/tree/master/vignettes/Intro.Rmd)
- [distIUPAC bam2IUPAC: fastq + reference >>> reference mapped bam + angsd >>> iupac fasta](https://github.com/kullrich/distIUPAC/tree/master/vignettes/bam2IUPAC.Rmd)
- [distIUPAC population distanceStats per sliding-window: iupac fasta + pop info >>> dStats ](https://github.com/kullrich/distIUPAC/tree/master/vignettes/dStats.Rmd)
- [distIUPAC genrate twisst input: iupac fasta >>> twisst trees](https://github.com/kullrich/distIUPAC/tree/master/vignettes/twisstTrees.Rmd)
- [distIUPAC on GTF featues: iupac fasta + GTF >>> feature specific distances ](https://github.com/kullrich/distIUPAC/tree/master/vignettes/GTFdistances.Rmd)

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
