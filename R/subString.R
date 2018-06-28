#' @title subString
#' @name subString
#' @description This function gets a subsequence from a \code{DNAString}, \code{RNAString}, \code{AAString}, \code{BString}, \code{DNAStringSet}, \code{RNAStringSet}, \code{AAStringSet}, \code{BStringSet} object from the \code{Biostrings} package.
#' @import Biostrings
#' @param x \code{DNAStringSet}, \code{RNAString}, \code{AAString}, \code{BString}, \code{DNAStringSet}, \code{RNAStringSet}, \code{AAStringSet}, \code{BStringSet}
#' @param s start vector
#' @param e end vector
#' @examples
#' data("MySequences", package = "distIUPAC")
#' myStart<-c(1,3)
#' myEnd<-c(2,4)
#' mySubSequence <- subString(MySequences, myStart, myEnd)
#' @seealso \link[XVector]{subseq}
#' @export subString
#' @author Kristian K Ullrich
subString<-function(x,s,e){
  x.class <- class(x)[1]
  se.matrix <- cbind(s,e)
  if(x.class=="DNAString"){
    myDNAString <- x
    newDNAString <- DNAString(paste(sapply(apply(se.matrix,1,function(x) subseq(myDNAString,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
    names(newDNAString) <- names(x)
    return(newDNAString)
  }
  if(x.class=="RNAString"){
    myRNAString <- x
    newRNAString <- RNAString(paste(sapply(apply(se.matrix,1,function(x) subseq(myRNAString,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
    names(newRNAString) <- names(x)
    return(newRNAString)
  }
  if(x.class=="AAString"){
    myAAString <- x
    newAAString <- AAString(paste(sapply(apply(se.matrix,1,function(x) subseq(myAAString,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
    names(newAAString) <- names(x)
    return(newAAString)
  }
  if(x.class=="BString"){
    myBString <- x
    newBString <- BString(paste(sapply(apply(se.matrix,1,function(x) subseq(myBString,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
    names(newBString) <- names(x)
    return(newBString)
  }
  if(x.class=="DNAStringSet"){
    myDNAStringSet <- x
    if(length(myDNAStringSet) > 1){
      newDNAStringSet <- DNAStringSet(apply(sapply(apply(se.matrix,1,function(x) subseq(myDNAStringSet,x[1],x[2])),function(x) paste0(x)),1,function(x) paste(x,sep="",collapse="")))
      names(newDNAStringSet) <- names(x)
      return(newDNAStringSet)
    }
    if(length(myDNAStringSet) == 1){
      newDNAStringSet <- DNAStringSet(paste(sapply(apply(se.matrix,1,function(x) subseq(myDNAStringSet,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
      names(newDNAStringSet) <- names(x)
      return(newDNAStringSet)
    }
  }
  if(x.class=="RNAStringSet"){
    myRNAStringSet <- x
    if(length(myRNAStringSet) > 1){
      newRNAStringSet <- RNAStringSet(apply(sapply(apply(se.matrix,1,function(x) subseq(myRNAStringSet,x[1],x[2])),function(x) paste0(x)),1,function(x) paste(x,sep="",collapse="")))
      names(newRNAStringSet) <- names(x)
      return(newRNAStringSet)
    }
    if(length(myRNAStringSet) == 1){
      newRNAStringSet <- RNAStringSet(paste(sapply(apply(se.matrix,1,function(x) subseq(myRNAStringSet,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
      names(newRNAStringSet) <- names(x)
      return(newRNAStringSet)
    }
  }
  if(x.class=="AAStringSet"){
    myAAStringSet <- x
    if(length(myAAStringSet) > 1){
      newAAStringSet <- AAStringSet(apply(sapply(apply(se.matrix,1,function(x) subseq(myAAStringSet,x[1],x[2])),function(x) paste0(x)),1,function(x) paste(x,sep="",collapse="")))
      names(newAAStringSet) <- names(x)
      return(newAAStringSet)
    }
    if(length(myAAStringSet) == 1){
      newAAStringSet <- AAStringSet(paste(sapply(apply(se.matrix,1,function(x) subseq(myAAStringSet,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
      names(newAAStringSet) <- names(x)
      return(newAAStringSet)
      }
  }
  if(x.class=="BStringSet"){
    myBStringSet <- x
    if(length(myBStringSet) > 1){
      newBStringSet <- BStringSet(apply(sapply(apply(se.matrix,1,function(x) subseq(myBStringSet,x[1],x[2])),function(x) paste0(x)),1,function(x) paste(x,sep="",collapse="")))
      names(newBStringSet) <- names(x)
      return(newBStringSet)
    }
    if(length(myBStringSet) == 1){
      newBStringSet <- BStringSet(paste(sapply(apply(se.matrix,1,function(x) subseq(myBStringSet,x[1],x[2])),function(x) paste0(x)),sep="",collapse=""))
      names(newBStringSet) <- names(x)
      return(newBStringSet)
    }
  }
}
