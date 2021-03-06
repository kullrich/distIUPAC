#' @title swgen
#' @name swgen
#' @description This function constructs a sliding window matrix
#' @param wlen window size [default: 1000]
#' @param wjump window jump [default: 1000]
#' @param start.by start position [default: 1]
#' @param end.by end position [default: 1000]
#' @examples
#' sw <- swgen(100, 100, 1, 1000)
#' @seealso \link[distIUPAC]{posgen}
#' @export swgen
#' @author Kristian K Ullrich
swgen<-function(wlen=1000, wjump=100, start.by=1, end.by=1000){
    if(end.by<start.by){
        stop("end.by < start.by")
    }
    start.seq<-seq(start.by, end.by, by = wjump)
    if((start.by + wlen - 1)>end.by){
        end.seq<-seq(end.by, end.by, by = wjump)
    }
    if((start.by + wlen - 1)<=end.by){
        end.seq<-seq(start.by + wlen - 1, end.by, by = wjump)
        end.seq<-c(end.seq, rep(end.by, length(start.seq)-length(end.seq)))
    }
    start.end.matrix<-rbind(start.seq, end.seq)
    start.end.matrix<-rbind(start.end.matrix, apply(start.end.matrix, 2,
      function(x) {mean( c( x[1] - 1, x[2]) )} ) )
    start.end.matrix[2,]<-as.numeric(sprintf("%.0f",start.end.matrix[2,]))
    row.names(start.end.matrix)<-c("start", "end", "mid")
    return(start.end.matrix)
}