#' @title swgen
#' @name swgen
#' @description This function constructs a sliding window matrix
#' @param wlen window size
#' @param wjump window jump
#' @param start.by start position
#' @param end.by end position
#' @examples
#' sw <- swgen(100, 100, 1, 1000)
#' @seealso \link[distIUPAC]{posgen}
#' @export swgen
#' @author Kristian K Ullrich
swgen<-function(wlen=1000,wjump=100,start.by=1,end.by=1000){
  if(end.by<=start.by){
    stop("end.by <= start.by")
  }
  start.seq<-seq(start.by, end.by, by = wjump)
  end.seq<-seq(start.by + wlen - 1, end.by, by = wjump)
  end.seq<-c(end.seq, rep(end.by, length(start.seq)-length(end.seq)))
  start.end.matrix<-rbind(start.seq, end.seq)
  start.end.matrix<-rbind(start.end.matrix, apply(start.end.matrix, 2, function(x) {mean( c( x[1] - 1, x[2]) )} ) )
  start.end.matrix[2,]<-as.numeric(sprintf("%.0f",start.end.matrix[2,]))
  row.names(start.end.matrix)<-c("start", "end", "mid")
  return(start.end.matrix)
}