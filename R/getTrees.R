#' @title getTrees
#' @name getTrees
#' @description This function returns nieghbor-joining (\link[ape]{njs} or \link[ape]{bionjs})
#' trees based on \code{distIUPAC} based distances rooted to a specific individual.
#' @import Biostrings
#' @import ape
#' @import doMC
#' @import foreach
#' @importFrom stats as.dist sd
#' @importFrom utils combn read.table setTxtProgressBar txtProgressBar
#' @param dna \code{DNAStringSet}
#' @param x.pos population X positions [default: NULL]
#' take all sequences from \code{DNAStringSet}
#' @param r.pos root position, needs to be within x.pos [default: 1]
#' @param resolve.root adds a zero-length branch below the MRCA
#' of the ingroup \link[ape]{root} [default: TRUE]
#' @param wlen sliding windows length  [default: 25000]
#' @param wjump sliding windows jump [default: 25000]
#' @param start.by optional start position [default: 1]
#' @param end.by optional end position [default: NULL]
#' @param wtype sliding windows type to use \code{bp}, \code{biSites}
#' or \code{triSites}  [default: "bp"]
#' @param dist distance to use [default: "IUPAC"]
#' @param global.deletion a logical indicating whether to delete the sites
#' with missing data in a global [default: TRUE] or pairwise way [FALSE]
#' @param model tree model to use either \link[ape]{njs}
#' or \link[ape]{bionjs} [default: "bionjs"]
#' @param threads number of parallel threads [default: 1]
#' @param ncores number of parallel cores to process pairwise distance
#' calculation [default: 1] see \link[distIUPAC]{rcpp_distIUPAC}
#' @param x.name population X name [default: "x"]
#' @param r.name root name [default: "r"]
#' @param chr.name chromosome name [default: "chr"]
#' @param plot indicates if tree should be plotted [default: FALSE]
#' @param pB specifies if progress should be shown as a progress bar
#' [default: FALSE]
#' @examples
#' data("MySequences", package = "distIUPAC")
#' x.pos<-c(1,2,3,4,5,35,63,71,79,82,88,89,97,105,106,114,115,116,117)
#' r.pos<-116
#' trees<-getTrees(MySequences, x.pos = x.pos, r.pos = r.pos,
#' wlen = 10000, wjump = 10000, x.name = "all", r.name = "Rnor", threads = 1)
#' trees
#' trees<-getTrees(MySequences, x.pos = x.pos, r.pos = r.pos,
#' wlen = 10000, wjump = 10000, x.name = "all", r.name = "Rnor",
#' threads = 2, plot = TRUE)
#' @export getTrees
#' @author Kristian K Ullrich
getTrees<-function(dna, x.pos=NULL, r.pos=1, resolve.root=TRUE, wlen=25000,
  wjump=25000, start.by=1, end.by=NULL, wtype="bp",
  dist="IUPAC", global.deletion=TRUE, model="bionjs", threads=1, ncores=1,
  x.name="x", r.name="r", chr.name="chr", plot=FALSE, pB=FALSE){
    options(scipen=22)
    if(!model%in%c("njs","bionjs")){
        stop("model needs to be either njs or bionjs")
    }
    if(length(x.pos)<3){stop("needs at least 3 sequences to calculate tree")}
    if(is.null(end.by)){end.by<-unique(width(dna))}
    if(start.by>unique(width(dna))){
        stop("start.by needs to be equal or smaller than dna length")
    }
    if(end.by>unique(width(dna))){
        stop("end.by needs to be equal or smaller than dna length")
    }
    if(wtype=="bp"){
        tmp.sw<-swgen(wlen=wlen, wjump=wjump, start.by=start.by, end.by=end.by)
    }
    if(wtype=="biSites"){
        tmp.POS<-biSites(dna, threads=threads, pB=FALSE)
        tmp.sw<-posgen(tmp.POS$biPOS, wlen=wlen, start.by=start.by,
          end.by=end.by)
    }
    if(wtype=="triSites"){
        tmp.POS<-triSites(dna, threads=threads, pB=FALSE)
        tmp.sw<-posgen(tmp.POS$triPOS, wlen=wlen, start.by=start.by,
          end.by=end.by)
    }
    if(is.null(x.pos)){
        x.pos<-seq(1, length(dna))
    }
    if(!r.pos%in%x.pos){x.pos<-c(x.pos, r.pos)}
    dna_<-dna[x.pos]
    x.pos_<-seq(1, length(x.pos))
    r.pos_<-which(r.pos==x.pos)
    j<-NULL
    if(pB){
        pb<-txtProgressBar(min=0, max=ncol(tmp.sw), initial=0, style=3)
    }
    registerDoMC(threads)
    OUT<-foreach(j=1:dim(tmp.sw)[2], .combine=rbind) %dopar% {
        XNAME<-x.name
        RNAME<-r.name
        CHRNAME<-chr.name
        START<-NA
        END<-NA
        tree.x<-NA
        topo.x<-NA
        treelength.x<-NA
        dSites.x<-NA
        dNA.x<-NA
        comment.x<-NA
        OUT<-list(XNAME, RNAME, CHRNAME, START, END, dSites.x,
          dNA.x, tree.x, topo.x, treelength.x, comment.x)
        names(OUT)<-c("XNAME", "RNAME", "CHRNAME", "START", "END", "dSites.x",
          "dNA.x", "tree.x", "topo.x", "treelength.x", "comment.x")
        OUT$START<-tmp.sw[1, j][[1]]
        OUT$END<-tmp.sw[2, j][[1]]
        tmp.seq<-subseq(dna_, OUT$START, OUT$END)
        if(global.deletion){
            tmp.seq<-globalDeletion(tmp.seq)
        }
        if(dist=="IUPAC"){
            dIUPAC<-rcpp_distIUPAC(as.character(tmp.seq), ncores=ncores)
            OUT$dSites.x<-mean(as.dist(dIUPAC$sitesUsed), na.rm=TRUE)
            OUT$dNA.x<-length(which(is.na(as.dist(
              dIUPAC$distIUPAC))))/length(as.dist(dIUPAC$distIUPAC))
        }
        if(dist!="IUPAC"){
            dIUPAC<-setNames( list(
              dist.dna(as.DNAbin(DNAMultipleAlignment(tmp.seq)),
                model=dist, as.matrix=TRUE, pairwise.deletion=TRUE),
              pairwiseDeletion(as.character(tmp.seq))$sitesUsed
                ), c("distIUPAC", "sitesUsed") )
            OUT$dSites.x<-mean(as.dist(dIUPAC$sitesUsed), na.rm=TRUE)
            OUT$dNA.x<-length(which(is.na(as.dist(
              dIUPAC$distIUPAC))))/length(as.dist(dIUPAC$distIUPAC))
        }
        if(model=="njs"){
            tmp.seq.tree<-tryCatch(njs(as.dist(dIUPAC$distIUPAC)),
              warning = function(war){return(paste0("WARNING: ", war))},
              error = function(err){return(paste0("ERROR: ", err))})
        }
        if(model=="bionjs"){
            tmp.seq.tree<-tryCatch(bionjs(as.dist(dIUPAC$distIUPAC)),
              warning = function(war){return(paste0("WARNING: ", war))},
              error = function(err){return(paste0("ERROR: ", err))})
        }
        if(class(tmp.seq.tree)!="phylo"){
            OUT$comment.x<-gsub("\n", "", tmp.seq.tree)
        }
        if(class(tmp.seq.tree)=="phylo"){
            OUT$tree.x<-write.tree(root(tmp.seq.tree, r.pos_,
              resolve.root=resolve.root))
            if(plot){
                tryCatch(
                  plot(root(tmp.seq.tree, r.pos_, resolve.root=resolve.root),
                    main=paste0(OUT$XNAME, " ", OUT$RNAME, "\n", OUT$CHRNAME,
                    "-", OUT$START, "-", OUT$END)),
                  warning = function(war){return(paste0("WARNING: ", war))},
                  error = function(err){return(paste0("ERROR: ", err))})
            }
            OUT$treelength.x<-sum(tmp.seq.tree$edge.length)
            tmp.seq.tree.topo<-tmp.seq.tree
            tmp.seq.tree.topo$edge.length<-NULL
            OUT$topo.x<-write.tree(root(tmp.seq.tree.topo, r.pos_,
              resolve.root=resolve.root))
        }
        if(pB){
            setTxtProgressBar(pb, j)
        }
        OUT
    }
    if(pB){
        setTxtProgressBar(pb, dim(tmp.sw)[2])
        close(pb)
    }
    return(OUT)
}
