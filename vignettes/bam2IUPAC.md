---
title: "01. bam2IUPAC"
author: "Kristian K Ullrich"
date: "2019-02-13"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{01. bam2IUPAC}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

In this vignette Users learn how to create  __`IUPAC fasta format`__ ([IUPAC code](https://www.bioinformatics.org/sms/iupac.html)) files from __`reference mapped bam`__ files which can be used with [distIUPAC](https://github.com/kullrich/distIUPAC) to calculate __`IUPAC distances`__.

## Pre-requistes for obtaining __fasta__ format files from reference mapped __bam__ files

The following external tools needs to be installed to be able to obtain IUPAC __fasta__ format files:

1. __`Genome mapper`__ of your choice (e.g. [bwa mem](http://bio-bwa.sourceforge.net/) or [NextGenMap](https://github.com/Cibiv/NextGenMap))
2. __`Picard tools`__ for __`bam`__ file sorting and de-duplication ([picard tools](http://broadinstitute.github.io/picard/))
3. __`angsd`__ for __`IUPAC fasta format`__ retrieval ([angsd](https://github.com/ANGSD/angsd))

A manual how the pre-requisites needs to be installed is given at the end of this vignette.

## STEP 1: Mapping __`fastq`__ files to a __`reference`__

### Using __`BWA`__ for mapping:

The user needs to perform five steps:

1. __`fastq`__ files, preferentially __`QC`__ pre-processed (see [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for one of many trimming tools available)
2. reference __`fasta`__ file
3. __`build index`__ for the reference
4. __`map reads`__
5. __`sort reads`__

```
#FORWARD FASTQ FILE: 1.fq
#REVERSE FASTQ FILE: 2.fq
#REFERENCE: ref.fasta
#USED THREADS: 12
#SAMPLE NAME: IND1
#SAM OUTPUT FILE: ind1.sam
#build index
bwa index ref.fasta
#map reads
bwa mem -t 12 -M -R @RG\tID:IND1\tLB:lib1\tPL:ILLUMINA\tPU:unit1\tSM:IND1 -o ind1.sam ref.fasta 1.fq 2.fq
#sort sam file
java -jar picard.jar SortSam I=ind1.sam O=ind1.sorted.bam SO=coordinate
```

### Using __`NextGenMap`__ for mapping

The usere needs to perform four steps:

1. __`fastq`__ files, preferentially __`QC`__ pre-processed (see [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) for one of many trimming tools available)
2. reference __`fasta`__ file
3. __`map reads`__
4. __`sort reads`__

```
#FORWARD FASTQ FILE: 1.fq
#REVERSE FASTQ FILE: 2.fq
#REFERENCE: ref.fasta
#USED THREADS: 12
#SAMPLE NAME: IND1
#SAM OUTPUT FILE: ind1.bam
#map reads
ngm -1 1.fq -2 2.fq -r ref.fasta -o ind1.bam --no-unal --sensitive -t 12 --no-progress --rg-id IND1 --rg-sm IND1 --rg-lb lib1 --rg-pl illumina --rg-pu unit1 -b
#sort bam file
java -jar picard.jar SortSam I=ind1.bam O=ind1.sorted.bam SO=coordinate
```
## OPTIONAL STEP: Merge __`sorted bam`__ files

In case that multiple __`fastq`__ libraries exists or have been mapped without prior merging of the __`fastq`__ files, the User can merge __`reference mapped bam`__ files as follows:

```
java -jar picard.jar MergeSamFiles I=ind1_1.sorted.bam I=ind1_2.sorted.bam O=ind1.sorted.bam
```

## STEP 2: Remove duplicates from __`sorted bam`__ file

```
java -jar picard.jar RemoveDuplicates REMOVE_DUPLICATES=true I=ind1.sorted.bam O=ind1.sorted.nodup.bam M=ind1.sorted.bam.duplicate.metrics
```

## STEP 3: Create __`IUPAC fasta`__ file for the de-duplicated __`referenced mapped bam`__ file

Chromosomes should be processed separately to be able to merge different samples into __`chromosome alignments`__, which can be processed with __`distIUPAC`__ as follows:

```
#chromosomes should be processed separately to be able to easily merge different samples into one alignment to be processed with 'distIUPAC'
#chromosome 'chr1' will be processed here
samtools index ind1.sorted.nodup.bam
angsd -doFasta 4 -doCounts 1 -minQ 20 -minMapQ 30 -uniqueOnly -setMinDepth 5 -setMaxDepth 100 -iupacRatio 0.2 -i ind1.sorted.nodup.bam -out IND1.minQ20.minMapQ30.uniqueOnly.setMinDepth5.setMaxDepth100.chr1 -r chr1
```

## STEP 4: Repeat STEP 1 to STEP 3 for multiple individuals

## STEP 5: Merge samples into chromosome alignments

Assuming all __`IUPAC fasta`__ files are located in the same folder and same chromosome files have the same ending, individuals can be merged as follows:

```
#chromosome 'chr1' from different individuals will be processed here
#1. uncompress fasta files
for file in *.chr1.fa.gz;do gunzip $file;done
#2. rename fasta sequences according to file names
for file in *.chr1.fa;do sed -i 's/>chr1/>'"$file"'/g' $file;done
#3. merge fasta files
for file in *.chr1.fa;do cat $file >> chr1.fa;done
```

## Calculating __`IUPAC distances`__ for 'chr1' with __`distIUPAC`__

The following commands needs to be executed in __`R`__:

```
library(distIUPAC)
dna<-readDNAStringSet("chr1.fa")
chr1.dist<-distIUPAC(as.character(subseq(dna,1,10000)))
```

## Pre-requisites installation - unix based systems (no MAC OS X)

Short description of how to compile the external tools for a unix based system are given.
See the next section for MAC OS X installation.

### NextGenMap installation
```
#download latest release from NextGenMap
wget https://codeload.github.com/Cibiv/NextGenMap/tar.gz/NextGenMap-0.5.5.tar.gz
tar -xvf v0.5.5.tar.gz
rm v0.5.5.tar.gz
cd NextGenMap-0.5.5/
mkdir -p build/
cd build/
cmake ..
make
```

### Picard tools installation
```
#download latest release 'picard.jar' from broadinstitute
wget https://github.com/broadinstitute/picard/releases/download/2.17.4/picard.jar
```

### ANGSD installation

```
#download latest release ANGSD from https://github.com/ANGSD/angsd
git clone https://github.com/samtools/htslib.git;
git clone https://github.com/angsd/angsd.git;
cd htslib;make;
cd ../angsd;
make HTSSRC=../htslib
```

## Pre-requisites installation - MAC OS X systems

For a MAC OS X system there are additional pre-requisites that needs to be installed to be able to compile all necessary software.

### Xcode Command Line Tools
'Xcode' needs to be installed from App-Store

### Homebrew
'Homebrew' needs to be installed [https://brew.sh/index_de.html](https://brew.sh/index_de.html)

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)
```

### Git
'Git' needs to be installed [https://git-scm.com/download/mac](https://git-scm.com/download/mac)

1. download latest Git DMG Image
2. Open DMG Iamge and Install PKG file

### Autoconf via Homebrew
'Autoconf' needs to be installed [http://www.gnu.org/software/autoconf/autoconf.html](http://www.gnu.org/software/autoconf/autoconf.html)
```
brew install autoconf
```

### Cmake
'Cmake' needs to be installed [https://cmake.org/download/](https://cmake.org/download/)

1. download latest cmake DMG Image release [https://cmake.org/files/v3.10/cmake-3.10.2-Darwin-x86_64.dmg](https://cmake.org/files/v3.10/cmake-3.10.2-Darwin-x86_64.dmg)
2. Open DMG Image and copy to Applications
3. in a Terminal
```
sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install
```

### NextGenMap installation
```
#download latest release from NextGenMap
curl -L https://github.com/Cibiv/NextGenMap/archive/v0.5.5.tar.gz > v0.5.5.tar.gz
tar -xvf v0.5.5.tar.gz
rm v0.5.5.tar.gz
cd NextGenMap-0.5.5/
mkdir -p build/
cd build/
cmake ..
make
```

### Picard tools installation
```
#download latest release 'picard.jar' from broadinstitute
curl -L https://github.com/broadinstitute/picard/releases/download/2.17.4/picard.jar > picard.jar
```

### ANGSD installation

```
#download latest release ANGSD from https://github.com/ANGSD/angsd
git clone https://github.com/samtools/htslib.git;
git clone https://github.com/angsd/angsd.git;
cd htslib;/usr/local/Cellar/autoconf/2.69/bin/autoconf;/usr/local/Cellar/autoconf/2.69/bin/autoheader;./configure --disable-lzma;make;
cd ../angsd;
make HTSSRC=../htslib
```

## References

Korneliussen TS, Albrechtsen A and Nielsen R. __ANGSD: Analysis of Next Generation Sequencing Data.__ _BMC Bioinformatics_ (2014) _15_:356 [https://doi.org/10.1186/s12859-014-0356-4](https://doi.org/10.1186/s12859-014-0356-4)

Li, H. & Durbin, R. __Fast and accurate short read alignment with Burrows-Wheeler transform.__ _Bioinformatics_ (2009) _25_:1754 [https://doi.org/10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324)

__Picard tools__ [http://broadinstitute.github.io/picard](http://broadinstitute.github.io/picard)

Sedlazeck FJ, Rescheneder P, von Haeseler A. __NextGenMap: fast and accurate read mapping in highly polymorphic genomes.__ _Bioinformatics_ (2013) _21_:2790 [https://doi.org/10.1093/bioinformatics/btt468](https://doi.org/10.1093/bioinformatics/btt468)
