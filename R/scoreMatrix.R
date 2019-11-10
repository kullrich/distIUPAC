#' @title scoreMatrix
#' @name scoreMatrix
#' @description This function creates a \code{scoreMatrix} object to be used
#' with the \code{distIUPACmatrix} function. By default, the score matrix is
#' defined as literal distance obtained from \code{Chang et al. 2017}.
#' (see \url{https://link.springer.com/article/10.1007/s00335-017-9704-9})
#' @references Chang, P. L., Kopania, E., Keeble, S., Sarver, B. A., Larson,
#' E., Orth, A., ... & Dean, M. D. (2017). Whole exome sequencing of
#' wild-derived inbred strains of mice improves power to link phenotype and
#' genotype. \emph{Mammalian genome}, \bold{28(9-10)}, 416-425.
#' @seealso \link[distIUPAC]{distIUPACmatrix}, \link[ape]{dist.dna}
#' @examples
#' data("MySequences", package="distIUPAC")
#' MyScoreMatrix <- scoreMatrix()
#' MyScoreMatrix["A","R"] <- 10.0
#' distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
#' MyScoreMatrix["A","R"] <- 0.5
#' distIUPACmatrix(as.character(MySequences[1:10]), MyScoreMatrix)
#' @export scoreMatrix
#' @author Kristian K Ullrich
scoreMatrix<-function(){
    distances<-c(
      ## A C G T  
      ## R Y S W K M  
      ## B D H V   
      ## . - N X
      # A
      0.0, 1.0, 1.0, 1.0,
      0.5, 1.0, 1.0, 0.5, 1.0, 0.5,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # C
      1.0, 0.0, 1.0, 1.0,
      1.0, 0.5, 0.5, 1.0, 1.0, 0.5,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # G
      1.0, 1.0, 0.0, 1.0,
      0.5, 1.0, 0.5, 1.0, 0.5, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # T
      1.0, 1.0, 1.0, 0.0,
      1.0, 0.5, 1.0, 0.5, 0.5, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # R
      0.5, 1.0, 0.5, 1.0,
      0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # Y
      1.0, 0.5, 1.0, 0.5,
      1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # S
      1.0, 0.5, 0.5, 1.0,
      1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # W
      0.5, 1.0, 1.0, 0.5,
      1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # K
      1.0, 1.0, 0.5, 0.5,
      1.0, 1.0, 1.0, 1.0, 0.0, 1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # M
      0.5, 0.5, 1.0, 1.0,
      1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # B
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # D
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # H
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # V
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # .
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # -
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # N
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      # X
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0,
      -1.0, -1.0, -1.0, -1.0
    )
    scoreMatrix<-matrix(distances, ncol=18, nrow=18)
    colnames(scoreMatrix)<-rownames(scoreMatrix)<-c(
      "A", "C", "G", "T", 
      "R", "Y", "S", "W", "K", "M",
      "B", "D", "H", "V",
      ".", "-", "N", "X")
    return(scoreMatrix)
}
