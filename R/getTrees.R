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
#' @param x.pos population X positions
#' @param r.pos root position, needs to be within x.pos
#' @param wlen sliding windows length
#' @param wjump sliding windows jump
#' @param start.by optional start position
#' @param end.by optional end position
#' @param wtype sliding windows type to use \code{bp}, \code{biSites} or \code{triSites}
#' @param dist distance to use
#' @param global.deletion a logical indicating whether to delete the sites with missing data in a global way (default is to delete in a pairwise way)
#' @param model tree model to use either \link[ape]{njs} or \link[ape]{bionjs}
#' @param threads number of parallel threads
#' @param x.name population X name
#' @param r.name root name
#' @param chr.name chromosome name
#' @param plot indicates if tree should be plotted
#' @examples
#' data("MySequences", package = "distIUPAC")
#' x.pos<-c(1,2,3,4,5,35,63,71,79,82,88,89,97,105,106,114,115,116,117)
#' r.pos<-116
#' trees<-getTrees(MySequences, x.pos = x.pos, r.pos = r.pos,
#' wlen = 10000, wjump = 10000, x.name = "all", r.name = "Rnor", threads = 1)
#' trees
#' trees<-getTrees(MySequences, x.pos = x.pos, r.pos = r.pos,
#' wlen = 10000, wjump = 10000, x.name = "all", r.name = "Rnor", threads = 1, plot = TRUE)
#' @export getTrees
#' @author Kristian K Ullrich
getTrees<-function(dna, x.pos=NULL, r.pos=1, wlen=25000, wjump=25000, start.by=NULL, end.by=NULL, wtype="bp", dist="IUPAC", global.deletion=FALSE, model="bionjs", threads=1, x.name="x", r.name="r", chr.name="chr", plot=FALSE){
  options(scipen=22)
  if(is.null(start.by)){start.by<-1}
  if(is.null(end.by)){end.by<-unique(width(dna))}
  if(start.by>unique(width(dna))){stop("start.by needs to be equal or smaller than dna length")}
  if(end.by>unique(width(dna))){stop("end.by needs to be equal or smaller than dna length")}
  if(!r.pos%in%x.pos){x.pos<-c(x.pos,r.pos)}
  if(length(x.pos)<3){stop("needs at least 3 sequences to calculate tree")}
  if(!model%in%c("njs","bionjs")){stop("model needs to be either njs or bionjs")}
  if(is.null(x.pos)){
    dna_<-dna
    x.pos_<-seq(1,length(dna))
  }
  if(!is.null(x.pos)){
    dna_<-dna[x.pos]
    x.pos_<-seq(1,length(x.pos))
  }
  x.pos_<-seq(1,length(x.pos))
  r.pos_<-which(r.pos==x.pos)
  if(wtype=="bp"){
    tmp.sw<-swgen(wlen=wlen,wjump=wjump,start.by=start.by,end.by=end.by)
  }
  if(wtype=="biSites"){
    tmp.POS<-biSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  if(wtype=="triSites"){
    tmp.POS<-triSites(dna_,x.pos_,threads=threads,pB=FALSE)
    tmp.sw<-posgen(tmp.POS,wlen=wlen,start.by=start.by,end.by=end.by)
  }
  j<-NULL
  pb<-txtProgressBar(min=0,max=dim(tmp.sw)[2],initial=0,style=3)
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
    OUT<-list(XNAME,RNAME,CHRNAME,START,END,dSites.x,dNA.x,tree.x,topo.x,treelength.x,comment.x)
    names(OUT)<-c("XNAME","RNAME","CHRNAME","START","END","dSites.x","dNA.x","tree.x","topo.x","treelength.x","comment.x")
    OUT$START<-tmp.sw[1,j][[1]]
    OUT$END<-tmp.sw[2,j][[1]]
    tmp.seq<-subseq(dna_,OUT$START,OUT$END)
    if(global.deletion){
      tmp.seq<-globalDeletion(tmp.seq)
    }
    if(dist=="IUPAC"){
      tmp.seq.dist<-distIUPAC(as.character(tmp.seq))
      OUT$dSites.x<-mean(as.dist(tmp.seq.dist$sitesUsed),na.rm=TRUE)
      OUT$dNA.x<-length(which(is.na(as.dist(tmp.seq.dist$distIUPAC))))/length(as.dist(tmp.seq.dist$distIUPAC))
      if(model=="njs"){
        tmp.seq.tree<-tryCatch(njs(as.dist(tmp.seq.dist$distIUPAC)), warning = function(war){return(paste0("WARNING: ",war))},
        error = function(err){return(paste0("ERROR: ",err))})
      }
      if(model=="bionjs"){
        tmp.seq.tree<-tryCatch(bionjs(as.dist(tmp.seq.dist$distIUPAC)), warning = function(war){return(paste0("WARNING: ",war))},
        error = function(err){return(paste0("ERROR: ",err))})
      }
      if(class(tmp.seq.tree)!="phylo"){
        OUT$comment.x<-gsub("\n","",tmp.seq.tree)
      }
      if(class(tmp.seq.tree)=="phylo"){
        OUT$tree.x<-write.tree(root(tmp.seq.tree,r.pos_))
        if(plot){plot(root(tmp.seq.tree,r.pos_),main=paste0(OUT$XNAME," ",OUT$RNAME,"\n",OUT$CHRNAME,"-",OUT$START,"-",OUT$END))}
        OUT$treelength.x<-sum(tmp.seq.tree$edge.length)
        tmp.seq.tree.topo<-tmp.seq.tree
        tmp.seq.tree.topo$edge.length<-NULL
        OUT$topo.x<-write.tree(root(tmp.seq.tree.topo,r.pos_))
      }
    }
    if(dist!="IUPAC"){
      tmp.seq.dist<-dist.dna(as.DNAbin.DNAMultipleAlignment(tmp.seq),model=dist,as.matrix=TRUE,pairwise.deletion=TRUE)
      tmp.seq.sites<-pairwiseDeletion(as.character(tmp.seq))$sitesUsed
      OUT$dSites.x<-mean(as.dist(tmp.seq.sites),na.rm=TRUE)
      OUT$dNA.x<-length(which(is.na(as.dist(tmp.seq.dist))))/length(as.dist(tmp.seq.dist))
      if(model=="njs"){
        tmp.seq.tree<-tryCatch(njs(as.dist(tmp.seq.dist)), warning = function(war){return(paste0("WARNING:  ",war))},
        error = function(err){return(paste0("ERROR:  ",err))})
      }
      if(model=="bionjs"){
        tmp.seq.tree<-tryCatch(bionjs(as.dist(tmp.seq.dist)), warning = function(war){return(paste0("WARNING:  ",war))},
        error = function(err){return(paste0("ERROR:  ",err))})
      }
      if(class(tmp.seq.tree)!="phylo"){
        OUT$comment.x<-gsub("\n","",tmp.seq.tree)
      }
      if(class(tmp.seq.tree)=="phylo"){
        OUT$tree.x<-write.tree(root(tmp.seq.tree,r.pos_))
        if(plot){plot(root(tmp.seq.tree,r.pos_),main=paste0(OUT$XNAME," ",OUT$RNAME,"\n",OUT$CHRNAME,"-",OUT$START,"-",OUT$END))}
        OUT$treelength.x<-sum(tmp.seq.tree$edge.length)
        tmp.seq.tree.topo<-tmp.seq.tree
        tmp.seq.tree.topo$edge.length<-NULL
        OUT$topo.x<-write.tree(root(tmp.seq.tree.topo,r.pos_))
      }
    }
    setTxtProgressBar(pb,j)
    OUT
  }
  setTxtProgressBar(pb,dim(tmp.sw)[2])
  close(pb)
  return(OUT)
}
