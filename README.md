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

### Quick-guide

### License

### Bug reports

Please report any errors or requests regarding `distIUPAC` to Kristian Ullrich (ullrich@evolbio.mpg.de)


