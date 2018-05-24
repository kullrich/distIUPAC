#' @title posgen
#' @name posgen
#' @description This function constructs a sliding window matrix based on a SNP position vector
#' @param pos position vector
#' @param wlen number of positions per window
#' @param start.by start position
#' @param end.by end position
#' @examples
#' sw <- posgen(seq(1, 975, by = 5), 50, 1, 1000)
#' @seealso \link[distIUPAC]{swgen}, \link[distIUPAC]{biSites}, \link[distIUPAC]{triSites}
#' @export posgen
#' @author Kristian K Ullrich
posgen<-function(pos,wlen=50,start.by=1,end.by=1000){
  if(end.by<=start.by){
    stop("end.by <= start.by")
  }
  start.seq<-seq(1, length(pos), by = wlen)
  end.seq<-seq(1 + wlen - 1, length(pos), by = wlen)
  end.seq<-c(end.seq, rep(length(pos), length(start.seq)-length(end.seq)))
  start.seq<-pos[start.seq]
  start.seq[1]<-start.by
  end.seq<-pos[end.seq]
  end.seq[length(end.seq)]<-end.by
  start.end.matrix<-rbind(start.seq, end.seq)
  for(i in 1:(dim(start.end.matrix)[2]-1)){
    new.end<-floor(mean(c(start.end.matrix[2,i],start.end.matrix[1,i+1])))
    new.start<-ceiling(mean(c(start.end.matrix[2,i],start.end.matrix[1,i+1])))
    start.end.matrix[2,i]<-new.end
    start.end.matrix[1,i+1]<-new.start
  }
  start.end.matrix<-rbind(start.end.matrix, apply(start.end.matrix, 2, function(x) {mean( c( x[1] - 1, x[2]) )} ) )
  start.end.matrix[2,]<-as.numeric(sprintf("%.0f",start.end.matrix[2,]))
  row.names(start.end.matrix)<-c("start", "end", "mid")
  return(start.end.matrix)
}